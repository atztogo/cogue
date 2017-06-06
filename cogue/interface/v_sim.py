from cogue.crystal.utility import get_oriented_lattice

#
# V_sim ascii
#    
def write_v_sim(cell, filename=None):
    lat = get_oriented_lattice(cell.lattice)
    text  = "# cogue generated file\n"
    # text += "%15.9f%15.9f%15.9f\n" % tuple(get_lattice_parameters(lattice))
    # text += "%15.9f%15.9f%15.9f\n" % tuple(get_angles(lattice))
    text += "%15.9f%15.9f%15.9f\n" % (lat[0,0], lat[0,1], lat[1,1])
    text += "%15.9f%15.9f%15.9f\n" % (lat[0,2], lat[1,2], lat[2,2])
    text += "#keyword: reduced\n"
    # text += "#keyword: angdeg, reduced\n"
    for s, p in zip(cell.get_symbols(), cell.get_points().T):
        text += "%15.9f%15.9f%15.9f %2s\n" % (p[0], p[1], p[2], s)

    if filename:
        w = open(filename, 'w')
        w.write(text)
        w.close()
    else:
        return text
