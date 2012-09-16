def line_plot( m, n, pt, color=None ):
    mlab.plot3d( [pt[m][0], pt[n][0]],
                 [pt[m][1], pt[n][1]],
                 [pt[m][2], pt[n][2]],
                 tube_radius=0.015, opacity=1, color=color )

def arrow_plot( lattice, color ):
    l = np.array([ x/np.linalg.norm(x) for x in lattice ])
    mlab.quiver3d( [ 0, 0, 0 ],
                   [ 0, 0, 0 ],
                   [ 0, 0, 0 ],
                   [ l[0,0], l[1,0], l[2,0] ],
                   [ l[0,1], l[1,1], l[2,1] ],
                   [ l[0,2], l[1,2], l[2,2] ],
                   color=color,
                   line_width=3,
                   scale_factor=1 )

    for c, v in zip( ('a','b','c'), l*1.3 ):
        mlab.text3d(v[0]+0.15, v[1], v[2], c, color=color, scale=0.3 )

def plot_lattice( lattice, color=None ):
    lat = lattice
    origin = np.zeros(3)
    pt = [ origin,
           lat[0],
           lat[1],
           lat[2],
           lat[0]+lat[1],
           lat[1]+lat[2],
           lat[2]+lat[0],
           lat[0]+lat[1]+lat[2] ]

    pt = np.array( pt )
    
    line_plot( 0, 1, pt, color )
    line_plot( 0, 2, pt, color )
    line_plot( 0, 3, pt, color )
    line_plot( 4, 7, pt, color )
    line_plot( 5, 7, pt, color )
    line_plot( 6, 7, pt, color )
    line_plot( 1, 4, pt, color )
    line_plot( 2, 5, pt, color )
    line_plot( 3, 6, pt, color )
    line_plot( 1, 6, pt, color )
    line_plot( 2, 4, pt, color )
    line_plot( 3, 5, pt, color )

def plot_axes( lattice, color ):
    arrow_plot( lattice, color )

def plot_lattice_points( lattice, dim ):
    positions = []
    for i in range( dim[0]+1 ):
        for j in range( dim[1]+1 ):
            for k in range( dim[2]+1 ):
                positions.append( [i,j,k] )

    p = np.dot( positions, lattice )
    mlab.points3d(p[:,0], p[:,1], p[:,2], scale_factor=0.2, opacity=0.2, color=(0,0,0))              

def plot_lattice_points_plusminus( lattice, dim ):
    positions = []
    for i in range( -dim[0], dim[0]+1 ):
        for j in range( -dim[1], dim[1]+1 ):
            for k in range( -dim[2], dim[2]+1 ):
                positions.append( [i,j,k] )

    p = np.dot( positions, lattice )
    mlab.points3d(p[:,0], p[:,1], p[:,2], scale_factor=0.2, opacity=0.2, color=(0,0,0))                    

def plot_modulation( lattice, positions, modulation ):
    modu = np.dot( modulation, lattice )
    x = positions[:,0]
    y = positions[:,1]
    z = positions[:,2]
    u = modu[:,0]
    v = modu[:,1]
    w = modu[:,2]
    mlab.quiver3d( x, y, z, u, v, w,
                   color=(0,0,0),
                   line_width=3,
                   scale_factor=10*options.amp_modulation )
    plot_atoms( positions )

def plot_atoms( positions, shift=[0,0,0] ):
    positions += shift
    x = positions[:,0]
    y = positions[:,1]
    z = positions[:,2]
    mlab.points3d(x, y, z, resolution=16, scale_factor=0.4, color=(0,0.5,0))                    

def fracval( frac ):
    if frac.find('/') == -1:
        return float( frac )
    else:
        x = frac.split('/')
        return float( x[0] ) / float( x[1] )

def set_shift( cell, shift ):
    scaled_positions = cell.get_scaled_positions()
    scaled_positions += shift
    scaled_positions -= np.floor( scaled_positions )
    cell.set_scaled_positions( scaled_positions )
