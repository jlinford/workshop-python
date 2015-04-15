import numpy as np
from mpi4py import MPI

# Initialize MPI
COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()
NPE = COMM.Get_size()
    

def advec_diff(cell_size, 
               c2l, w2l, d2l, 
               c1l, w1l, d1l,
               c, w, d,
               c1r, w1r, d1r,
               c2r, w2r, d2r):
    """
    Upwinded advection/diffusion equation.
    c = conc, w = wind, d = diff
    x2l is the 2-left neighbor of x, etc.
    x2r is the 2-right neighbor of x, etc.
    """
    wind = (w1l + w) / 2.0
    if wind >= 0.0:
        advec_termL = (1.0/6.0) * ( -c2l + 5.0*c1l + 2.0*c )
    else:
        advec_termL = (1.0/6.0) * ( 2.0*c1l + 5.0*c - c1r )
    advec_termL *= wind
    wind = (w1r + w) / 2.0
    if wind >= 0.0: 
        advec_termR = (1.0/6.0) * ( -c1l + 5.0*c + 2.0*c1r )
    else:
        advec_termR = (1.0/6.0) * ( 2.0*c + 5.0*c1r - c2r )
    advec_termR *= wind
    advec_term = (advec_termL - advec_termR) / cell_size
    diff_term = ( ((d1l+d)/2.0)*(c1l-c) - ((d+d1r)/2.0)*(c-c1r) ) / (cell_size * cell_size)
    return advec_term + diff_term


def space_advec_diff(cell_size, c, w, d, cb, wb, db, dcdx):
    """
    Apply advection / diffusion
    """
    n = len(c)
    
    # Boundary cell c[0]
    dcdx[0] = advec_diff(cell_size,
                         cb[0], wb[0], db[0],  # 2-left neighbors
                         cb[1], wb[1], db[1],  # 1-left neighbors
                         c[0], w[0], d[0],     # Values
                         c[1], w[1], d[1],     # 1-right neighbors
                         c[2], w[2], d[2])     # 2-right neighbors

    # Boundary cell c[1]
    dcdx[1] = advec_diff(cell_size,
                         cb[1], wb[1], db[1],  # 2-left neighbors
                         cb[2], wb[2], db[2],  # 1-left neighbors
                         c[1], w[1], d[1],     # Values
                         c[2], w[2], d[2],     # 1-right neighbors
                         c[3], w[3], d[3])     # 2-right neighbors

    for i in xrange(2,n-2):
        dcdx[i] = advec_diff(cell_size,
                             c[i-2], w[i-2], d[i-2],  # 2-left neighbors
                             c[i-1], w[i-1], d[i-1],  # 1-left neighbors
                             c[i],   w[i],   d[i],    # Values
                             c[i+1], w[i+1], d[i+1],  # 1-right neighbors
                             c[i+2], w[i+2], d[i+2])  # 2-right neighbors

    # Boundary cell c[n-2]
    dcdx[n-2] = advec_diff(cell_size,
                           c[n-4], w[n-4], d[n-4],  # 2-left neighbors
                           c[n-3], w[n-3], d[n-3],  # 1-left neighbors
                           c[n-2], w[n-2], d[n-2],  # Values
                           cb[1],  wb[1],  db[1],   # 1-right neighbors
                           cb[2],  wb[2],  db[2])   # 2-right neighbors

    # Boundary cell c[n-1]
    dcdx[n-1] = advec_diff(cell_size,
                           c[n-3], w[n-3], d[n-3],  # 2-left neighbors
                           c[n-2], w[n-2], d[n-2],  # 1-left neighbors
                           c[n-1], w[n-1], d[n-1],  # Values
                           cb[2],  wb[2],  db[2],   # 1-right neighbors
                           cb[3],  wb[3],  db[3])   # 2-right neighbors


def discretize(cell_size, dt, c0, w, d, cb, wb, db):
    """
    Discretization kernel
    """
    out = np.copy(c0)
    c = np.copy(c0)
    dcdx = np.empty_like(c0)

    space_advec_diff(cell_size, c0, w, d, cb, wb, db, dcdx)
    c = c + dt * dcdx
    space_advec_diff(cell_size, c, w, d, cb, wb, db, dcdx)
    c = c + dt * dcdx
    
    return map(lambda x: max(x, 0), 0.5 * (out + c))


def discretize_rows(dx, dt, conc, wind, diff):
    cb = np.empty(4)
    wb = np.empty(4)
    db = np.empty(4)
    m, n = conc.shape
    for i in xrange(m):
        cb[0] = conc[i][n-2]
        cb[1] = conc[i][n-1]
        cb[2] = conc[i][0]
        cb[3] = conc[i][1]
        wb[0] = wind[i][n-2]
        wb[1] = wind[i][n-1]
        wb[2] = wind[i][0]
        wb[3] = wind[i][1]
        db[0] = diff[i][n-2]
        db[1] = diff[i][n-1]
        db[2] = diff[i][0]
        db[3] = diff[i][1]
        conc[i,:] = discretize(dx, dt, conc[i,:], wind[i,:], diff[i,:], cb, wb, db)


def discretize_cols(dy, dt, conc, wind, diff):
    cb = np.empty(4)
    wb = np.empty(4)
    db = np.empty(4)
    m, n = conc.shape
    
    shape = (2, n)
    north_conc = np.empty(shape)
    south_conc = np.empty(shape)
    north_wind = np.empty(shape)
    south_wind = np.empty(shape)
    north_diff = np.empty(shape)
    south_diff = np.empty(shape)

    north = NPE - 1 if RANK == 0 else RANK - 1
    south = 0 if RANK == (NPE-1) else RANK + 1

    COMM.Sendrecv(sendbuf=conc[:2,:], dest=north, sendtag=11,
                  recvbuf=south_conc, source=south, recvtag=11)
    COMM.Sendrecv(sendbuf=conc[-2:,:], dest=south, sendtag=12,
                  recvbuf=north_conc, source=north, recvtag=12)
    COMM.Sendrecv(sendbuf=wind[:2,:], dest=north, sendtag=11,
                  recvbuf=south_wind, source=south, recvtag=11)
    COMM.Sendrecv(sendbuf=wind[-2:,:], dest=south, sendtag=12,
                  recvbuf=north_wind, source=north, recvtag=12)
    COMM.Sendrecv(sendbuf=diff[:2,:], dest=north, sendtag=11,
                  recvbuf=south_diff, source=south, recvtag=11)
    COMM.Sendrecv(sendbuf=diff[-2:,:], dest=south, sendtag=12,
                  recvbuf=north_diff, source=north, recvtag=12)

    for j in xrange(n):
        cb[0] = north_conc[0][j]
        cb[1] = north_conc[1][j]
        cb[2] = south_conc[0][j]
        cb[3] = south_conc[1][j]
        wb[0] = north_wind[0][j]
        wb[1] = north_wind[1][j]
        wb[2] = south_wind[0][j]
        wb[3] = south_wind[1][j]
        db[0] = north_diff[0][j]
        db[1] = north_diff[1][j]
        db[2] = south_diff[0][j]
        db[3] = south_diff[1][j]
        conc[:,j] = discretize(dy, dt, conc[:,j], wind[:,j], diff[:,j], cb, wb, db)


def fixedgrid(nrows=50, ncols=50, 
              source_rate=4.67e21, source_row=25, source_col=25,
              dx=1000, dy=1000, dt=50,
              tstart=0, tend=3600,
              O3_init=8.61e09,
              wind_init=(5.0, -15.0),
              diff_init=100):
    """
    Simple upwind finite differences fixed grid advection/diffusion model with dimension split X/Y
    """
    my_cols = ncols
    my_rows = nrows / NPE
    if RANK == (NPE-1):
        my_rows += (nrows % NPE)
    print '%d: %d rows' % (RANK, my_rows)
    my_row0 = (nrows / NPE) * RANK
    
    shape = (my_rows, my_cols)
    source = (source_row, source_col)
    conc = np.full(shape, O3_init)
    wind_u = np.full(shape, wind_init[0])
    wind_v = np.full(shape, wind_init[1])
    diff = np.full(shape, diff_init)
    time = tstart
    
    if (source_row > my_row0) and (source_row < (my_row0+my_rows)):
        print 'Rank %d added plume at %s' % (RANK, str(source))
        conc[source_row-my_row0,source_col] += source_rate / (dx * dy)

    if RANK == 0:
        print 'TIME: %g' % time
    while time < tend:
        discretize_rows(dx, dt/2.0, conc, wind_u, diff)
        discretize_cols(dy, dt, conc, wind_v, diff)
        discretize_rows(dx, dt/2.0, conc, wind_u, diff)
        time += dt
        if RANK == 0:
            print 'TIME: %g' % time
    if RANK == 0:
        print 'done'


if __name__ == "__main__":
  fixedgrid(tstart=0, tend=900)



