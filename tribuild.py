'''
      
                Various topology manipulation tools that can handle arbitrary simulation box shapes

                      Anze Hubman, Theory dept./National Institute of Chemistry, Slovenia
             
                                               August 2022

'''

'''
    This library is intended to handle general (e.g. triclinic) simulation box shapes which are
    usually not supported in standard topology builders.

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

# def find_angles(type_i,type_j,type_k,atom_id,element,x,y,z,r_bond_ij,r_bond_jk):
#    ''' declares an angle ijk (j is the central atom) if d(i,j) <= r_bond_ij and d(j,k) <= r_bond_jk '''


                            
                            
    
    

        
    

    
    

    

    

