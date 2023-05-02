from fipy.meshes import Grid3D
import numpy as np
import sys
import os


ndiv = int(sys.argv[1])
directory = sys.argv[2]
if not os.path.isdir(directory):
    os.mkdir(directory)

print('ndiv',ndiv)
nx = ndiv
ny = ndiv
nz = ndiv




baseMesh = Grid3D(dx = 1.0 / nx , dy = 1.0 / ny, dz = 1.0 / nz,
                  nx = nx, ny = ny, nz = nz)

coord=baseMesh.vertexCoords
np.savetxt(directory+'/coord.dat', coord.transpose(),
           fmt='%.10e %.10e %.10e', delimiter=' ', newline='\n',
           header='', footer='', comments='# ')


#print('facevertex')
facevertex=baseMesh.faceVertexIDs
#print(facevertex)
#print(' ' )
np.savetxt(directory+'/faces.dat', facevertex.transpose()+1,
           fmt='%d %d %d %d', delimiter=' ', newline='\n',
           header='', footer='', comments='# ')

#print('cellfaces')
cellFaceIDs = baseMesh.cellFaceIDs
#print(cellFaceIDs)
#print(' ')
ncell=cellFaceIDs.shape[1]
np.savetxt(directory+'/faces_in_cell.dat', cellFaceIDs.transpose()+1,
           fmt='%d %d %d %d %d %d', delimiter=' ', newline='\n',
           header='', footer='', comments='# ')

topol=np.zeros([8,ncell],dtype=np.int32)
for i in range(cellFaceIDs.shape[1]):
    c = cellFaceIDs[:,i]
    #print(i,' ',c)
    nodes = facevertex[:,c]
    vertices = np.unique(nodes.flatten())
    #print('nodes_in_faces')
    #print(nodes)
    #print('nodes_in_cell')
    #print(vertices)
    topol[:,i]=vertices
#print(' ')

np.savetxt(directory+'/topol.dat', topol.transpose()+1,
           fmt='%d %d %d %d %d %d %d %d', delimiter=' ', newline='\n',
           header='', footer='', comments='# ')
    
    
exterior,  =np.where(baseMesh.exteriorFaces)
#print(exterior)

interior, =np.where(baseMesh.interiorFaces)
#print(interior)

np.savetxt(directory+'/internal.dat', interior+1,
           fmt='%d', delimiter=' ', newline='\n',
           header='', footer='', comments='# ')
