'''
      
                Various topology manipulation tools that can handle arbitrary simulation box shapes

                      Anze Hubman, Theory dept./National Institute of Chemistry, Slovenia
             
                                               August 2022

'''

'''
    This module is intended to handle general (e.g. triclinic) simulation box shapes.

    Important note: this is currently development version. Several routines
                    will be added and extensively tested.
                    To contribute or report a bug please me an e-mail: anze.hubman@ki.si.
                    Use and modify freely but with caution. I do not take responsibility
                    for any mistakes that may occur.

'''  



import numpy as np

def generate_H_matrix(ax,ay,az,bx,by,bz,cx,cy,cz):
    ''' generate H matrix - parallelepiped for arbitrary sim. boxes '''
    
    H = np.zeros((3,3))
                           # General parallelepiped:
    H[0][0] = ax           # [ ax bx cx ]
    H[0][1] = bx           # [ ay by cy ]
    H[0][2] = cx           # [ az bz cz ]
    H[1][0] = ay           # parameters can be easily obtained from .cfg file or from
    H[1][1] = by           # a .cif file using e.g. Atomsk
    H[1][2] = cy
    H[2][0] = az
    H[2][1] = bz
    H[2][2] = cz

    return H

def invert_H_matrix(H):
    ''' compute inverse of H matrix '''
    
    inv_H = np.linalg.inv(H)
    return inv_H

def volume(H):
    ''' compute volume of a parallelepiped '''

    V = np.det(H)
    return V

def distance_PBC(atom_i,atom_j,H):
    ''' compute distance under PBC for arbitrary sim. boxes '''

    inv_H = invert_H_matrix(H)

    # convert to fractional coordinates
    atom_i_frac = np.dot(inv_H,atom_i)
    atom_j_frac = np.dot(inv_H,atom_j)
    s_ij_frac   = atom_i_frac - atom_j_frac

    # apply minimum image convention
    s_ij_frac[0] = s_ij_frac[0] - np.rint(s_ij_frac[0])
    s_ij_frac[1] = s_ij_frac[1] - np.rint(s_ij_frac[1])
    s_ij_frac[2] = s_ij_frac[2] - np.rint(s_ij_frac[2])

    # convert fractional to Cartesian
    s_ij_cart = np.dot(H,s_ij_frac)

    # compute distance
    R = np.linalg.norm(s_ij_cart)

    return R

def read_xyz(filename_xyz):
    ''' read atomic coordinates from .xyz file '''

    atom_id = []     # string
    element = []     # string
    x       = []     # float
    y       = []     # float
    z       = []     # float

    f = open(filename_xyz,'r')

    line_number = 0
    
    for line in f:
        
        if (line_number > 1):
            q = line.split(' ')
            l = []

            for elt in q:
                if (elt != '') and (elt != '\n'):
                    l.append(elt)

            atom_id.append(l[0])
            element.append(l[1])
            x.append(float(l[2]))
            y.append(float(l[3]))
            z.append(float(l[4]))

        line_number += 1

    f.close()

    return atom_id, elemnt, x, y, z


def find_bonds(type_i,type_j,atom_id,element,x,y,z,r_bond):
    ''' declares a bond between two atoms of types type_i and type_j if d(i,j) <= r_bond '''

    n = len(atom_id)
    bonds = []

    for i in range(n):

        atom_i = [x[i],y[i],z[i]]

        if (element[i] == type_i):

            for j in range(n):

                atom_j = [x[j],y[j],z[j]]

                if (element[j] == type_j):

                    r = distance_pbc(atom_i,atom_j,H)

                    if (type_i != type_j):

                        if (r <= r_bond):
                            tmp = np.sort([atom_id[i],atom_id[j]])
                            bonds.append(tmp.tolist())

                    if (type_i == type_j):

                        if ((r <= r_bond) and (i != j)):
                            tmp = np.sort([atom_id[i],atom_id[j]])      # sorting ensures there are no duplicates
                            tmp = tmp.tolist()                          # when type_i = type_j

                            if (tmp not in bonds):
                                bonds.append(tmp)

    return bonds

def read_poscar(filename='POSCAR',seldyn='no'):
    ''' reads VASP POSCAR file, returns H matrix and Cartesian atomic coordinates;
        seldyn declares whether Selective dynamics tag is present or not '''

    # read POSCAR file
    f = open(filename,'r')
    tmp = []
    for line in f:
        q = line.strip('\n').split(' ')
        l = []
        for elt in q:
            if (elt != '') and (elt != '\n'):
                l.append(elt)
        tmp.append(l)

    # scaling factor
    s = float(tmp[1][0])

    # generate H matrix
    ax = float(tmp[2][0])
    ay = float(tmp[2][1])
    az = float(tmp[2][2])
    bx = float(tmp[3][0])
    by = float(tmp[3][1])
    bz = float(tmp[3][2])
    cx = float(tmp[4][0])
    cy = float(tmp[4][1])
    cz = float(tmp[4][2])

    H = generate_H_matrix(ax,ay,az,bx,by,bz,cx,cy,cz)

    # species
    species_type = tmp[5]
    species_nmbr = tmp[6]

    # coordinates
    if (seldyn == 'no'):
        crd_type = tmp[7][0]
        i_start  = 8
    if (seldyn == 'yes'):
        crd_type = tmp[8][0]
        i_start  = 9
    
    x = []
    y = []
    z = []

    for i in range(i_start,len(tmp)):
        x_tmp = float(tmp[i][0])
        y_tmp = float(tmp[i][1])
        z_tmp = float(tmp[i][2])

        if (crd_type == 'Cartesian'):
            x.append(s*x_tmp)
            y.append(s*y_tmp)
            z.append(s*z_tmp)

        if (crd_type == 'Direct'):
            frac = [x_tmp,y_tmp,z_tmp]
            cart = np.dot(H,frac)
            x.append(cart[0])
            y.append(cart[1])
            z.append(cart[2])

    return x,y,z,H

def radial_shift(x,y,z,r,H,fixed_index,move_index_start,move_index_end,move_index):
    ''' shifts a moiety given in the topology file along a specified vector '''

    # unit vector from fixed atom to chosen atom of a moving moiety 
    pos_move  = []
    for i in range(move_index_start-1,move_index_end):
        pos_move.append([x[i],y[i],z[i]])

    pos_move  = np.array(pos_move)
    pos_fixed = np.array([x[fixed_index-1],y[fixed_index-1],z[fixed_index-1]])
    ori_move  = np.array([x[move_index-1], y[move_index-1],z[move_index-1]])

    r_0 = np.linalg.norm(ori_move-pos_fixed)

    e = (ori_move - pos_fixed)/r_0

    # shift molecule in the radial direction
    shift_x = e[0]*(r-r_0)
    shift_y = e[1]*(r-r_0)
    shift_z = e[2]*(r-r_0)

    n_move = move_index_end - move_index_start + 1
    for i in range(n_move):
        pos_move[i][0] += shift_x
        pos_move[i][1] += shift_y
        pos_move[i][2] += shift_z

    # checks:
    # a) if the final distance is correct
    # b) if PBC affects the minimum distance. Note that currently shifting does not obey
    #    the minimimum image convention. A more general implementation is planned. 
    r_check = np.linalg.norm(pos_move[0] - pos_fixed)
    r_check_pbc = distance_PBC(pos_move[0],pos_fixed,H)

    print(r_0,r_check,r_check_pbc,np.abs(r_check-r_check_pbc))

    return pos_move


                            
    
    

        
    

    
    

    

    

