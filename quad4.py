#!/usr/bin/env

# Translation of quad4.m from Matlab to Python
# by Carl Sandrock carl.sandrock@gmail.com
# 2014-03-15

from __future__ import division
from matplotlib.pylab import *
from scipy.sparse import coo_matrix
import time

def readmatrix(f, N, types, rows=False):
    returns = []
    for i in range(N):
        items = f.next().split()
        returns.append([convert(item) for convert, item in zip(types, items)])
    if rows:
        return returns
    else:
        return zip(*returns)

filename = 'sample.inp'
file_out = 'sample.out'

fid = open(filename, 'r')

tic = time.time()
# Read number of nodes
fid.next()
nnodes = int(fid.next())

#Read nodal coordinates
fid.next()
nodes, c1, c2 = readmatrix(fid, nnodes, [int, float, float])
coor = array(zip(c1, c2))

# Read number of elements
fid.next()
nelem = int(fid.next())

# Read if elements are plane stress or plain strain
fid.next()
plane = bool(fid.next())

# Read element number and element connectivity 
fid.next()
elnodes = array(readmatrix(fid, nelem, [int]*5, rows=True))

# Read material constants and element thickness
fid.next()
[[elas], [pois], [t]] = readmatrix(fid, 1, [float]*3)

# Read number of prescribed displacements
fid.next()
ndispl = int(fid.next())

# Read prescribed displacements
fid.next()
displ = array(sorted(readmatrix(fid, ndispl, [int, float, float], rows=True)))

# Read number of nodal loads
fid.next()
ncload = int(fid.next())

# Read nodal loads
fid.next()
cload = readmatrix(fid, ncload, [float]*3, rows=True)
fid.close()

finish = time.time() - tic
print 'Done reading input file         : ', finish, ' seconds'

tic = time.time()

# Find specified (sdof) and unknown (udof) degrees of freedom
dof = ones([nnodes, 2])
Ub = matrix(displ[:,2]).T
for node, dimension, displacement in displ:
	dof[nodes.index(node), dimension-1] = 0
sdof = flatnonzero(dof == 0)
udof = flatnonzero(dof != 0)

# Construct D matrix
# If plane strain
if plane == 0:
    e = elas/(1 - pois**2)
    nu = pois/(1 - pois)
else:
    e = elas
    nu = pois

# Not yet checked
c = e/(1 - nu**2)
E = matrix([[c,    c*nu,            0],
	        [c*nu, c,               0],
	        [0,    0,  0.5*c*(1 - nu)]])

# #Initialize global stiffness matrix
row_vec   = zeros(64*nelem)
col_vec   = zeros(64*nelem)
stiff_vec = zeros(64*nelem)
 
# #Initialize global load vector
p_global = zeros(2*nnodes)

# #Initialize B_ALL and DB_ALL
B_ALL  = zeros((12*nelem, 8))
EB_ALL = zeros((12*nelem, 8))

Gauss_pos = 1.0/sqrt(3.0)
# #Main loop over elements. Compute k_elem and assemble
B     = asmatrix(zeros((3, 8)))
D     = asmatrix(zeros((4, 8)))
dNdxi = asmatrix(zeros((2, 4)))
B2    = asmatrix(zeros((3, 4)))
invJ  = asmatrix(zeros((2, 2)))
pg    = zeros(8)

pos_vec = 0
for i in range(nelem):
    # Find coordinates of element nodes
    x = coor[elnodes[i, 1:5] - 1, 0]
    y = coor[elnodes[i, 1:5] - 1, 1]
    XY = matrix(zip(x, y))
    # Numerical integration loops    
    k_elem  = zeros((8,8))

    Gausspt = 0
    for jGauss in range(2):
        for iGauss in range(2):
            xi  = (-1)**(iGauss + 1)*Gauss_pos
            eta = (-1)**(jGauss + 1)*Gauss_pos
            dNdxi[:] = [[0.25*(eta - 1.0), 0.25*(1.0 - eta), 0.25*(1.0 + eta), 0.25*(-1.0 - eta)],
                        [0.25*(xi - 1.0), 0.25*(-1.0 - xi), 0.25*(1.0 + xi), 0.25*(1.0 - xi)]]
            D[:2,::2] = dNdxi
            D[2:,1::2] = dNdxi
            J = dNdxi*XY
            detJ = det(J)
            invJ = J.I
            B2[1 - 1, 1 - 1] = invJ[1 - 1, 1 - 1]
            B2[1 - 1, 2 - 1] = invJ[1 - 1, 2 - 1]
            B2[2 - 1, 3 - 1] = invJ[2 - 1, 1 - 1]
            B2[2 - 1, 4 - 1] = invJ[2 - 1, 2 - 1]
            B2[3 - 1, 1 - 1] = invJ[2 - 1, 1 - 1]
            B2[3 - 1, 2 - 1] = invJ[2 - 1, 2 - 1]
            B2[3 - 1, 3 - 1] = invJ[1 - 1, 1 - 1]
            B2[3 - 1, 4 - 1] = invJ[1 - 1, 2 - 1]
            B = B2*D
            k_elem += t*B.T*E*B*detJ
            B_ALL[i*12 + Gausspt*3:i*12+Gausspt*3 + 3, :] = B
            EB_ALL[i*12 + Gausspt*3:i*12+Gausspt*3 + 3, :] = E*B
            Gausspt += 1

    # Assemble k_elem into k_global
    p     = elnodes[i,1:5]
    pg[:] = [2*p[1 - 1]-1,
             2*p[1 - 1],
             2*p[2 - 1]-1,
             2*p[2 - 1],
             2*p[3 - 1]-1,
             2*p[3 - 1],
             2*p[4 - 1]-1,
             2*p[4 - 1]]
    
    for ii in range(8):
        for jj in range(8):
            if True or pg[jj] >= pg[ii]:
                row_vec[pos_vec]   = pg[ii]-1
                col_vec[pos_vec]   = pg[jj]-1
                stiff_vec[pos_vec] = k_elem[ii, jj]
                pos_vec += 1

# end of main loop over elements
k_global = coo_matrix((stiff_vec, (row_vec, col_vec)),
                      shape=(2*nnodes, 2*nnodes))
# FIXME: k_global and the other matrices based on it is not sparse!
k_global = k_global.todense().T

# clear pos_vec row_vec stiff_vec

# Add nodal loads to global load vector
for node, dimension, p_load in cload:
    p_global[nodes.index(node)*2 + dimension - 1] += p_load

# #Partition k_global into Kaa, Kab, Kba, Kbb
Kaa = k_global[np.ix_(udof, udof)]
Kab = k_global[np.ix_(udof, sdof)]
Kbb = k_global[np.ix_(sdof, sdof)]
#clear k_global

# # Modify RHS to include specified displacements
Pa  = p_global[udof, None] - Kab*Ub

finish = time.time() - tic
print 'Done assembling stiffness matrix:', finish, 'seconds.'

tic = time.time()
#Solve unknown dof's
# Note this does the lower factorisation while Matlab/Octave does upper by default
LKaa = linalg.cholesky(Kaa)
temp = solve(LKaa, Pa)
Ua = solve(LKaa.T, temp)

finish = time.time() - tic
print 'Done solving system:', finish, 'seconds.'

tic = time.time()

# Solve support reactions
Pb = Kab.T*Ua + Kbb*Ub

# # Sort Ua and Ub into A
U = matrix(zeros((2*nnodes, 1)))
U[udof, :] = Ua
U[sdof, :] = Ub

strain = zeros((nelem, nnodes))
stress = zeros((nelem, nnodes))

#Compute element strains and stresses
for i in range(nelem):
    p  = elnodes[i,1:5]
    pg = 2*array([p[0]-0.5, p[0], p[1]-0.5, p[1], p[2]-0.5, p[2], p[3]-0.5, p[3]], int) - 1
    strain[i, :] = (B_ALL[i*12:12*(i+1), :].dot(U[pg, :])).T
    stress[i, :] = (EB_ALL[i*12:12*(i+1), :].dot(U[pg, :])).T

# Write output to file
# Uoutput = hstack([array(nodes)[:, None], U[::2, :], U[1::2, :]])
# StressOut = hstack([elnodes[:, 1, None], stress[:, 0:6], stress[:, 9:12], stress[:, 6:9]])

# StressNode = zeros((nelem, 13))
# StressNode[:,0] = elnodes[:, 0]

# factors = array([1 + sqrt(3)/2, -0.5, 1., (1 - sqrt(3))/2])
# columns = array([[1,  4, 10, 7],
#                  [4,  1,  7, 10],
#                  [7,  4, 10, 1],
#                  [10, 1,  7, 4]])
# for i in range(3):
#     for column in columns:
#         print [i + c + 1 for c in column]
#         StressNode[:, i + column[0]] = (factors * StressOut[:, column+i]).sum(1)

# VonMises = zeros(nelem,4)
# if (plane==1)
#     for i=1:4
#         VonMises(:,i) = StressNode(:,2+(i-1)*3).**2 - ...
#                 StressNode(:,2+(i-1)*3).*StressNode(:,3+(i-1)*3) + ...
#                 StressNode(:,3+(i-1)*3).**2 + 3*StressNode(:,4+(i-1)*3).**2
#         VonMises(:,i) = VonMises(:,i).**0.5
#     end
# else
#     for i=1:4
#         VonMises(:,i) = (1-pois+pois**2)*(StressNode(:,2+(i-1)*3).**2 + ...
#                         StressNode(:,3+(i-1)*3).**2) - (1+pois-pois**2)* ...
#                         StressNode(:,2+(i-1)*3).*StressNode(:,3+(i-1)*3)+...
#                         3*StressNode(:,4+(i-1)*3).**2
#         VonMises(:,i) = VonMises(:,i).**0.5
#     end
# end    

# Tresca = zeros(nelem,4)
# if plane ==1
#     for j=1:nelem
#         for i=1:4
#             s = [StressNode(j,2+(i-1)*3) StressNode(j,4+(i-1)*3) 0
#                  StressNode(j,4+(i-1)*3) StressNode(j,3+(i-1)*3) 0
#                 0 0 0]
#             principal = eig(s)
#             Tresca(j,i) = max(principal)-min(principal)
#         end
#     end
# else
#     for j=1:nelem
#         for i=1:4
#             s = [StressNode(j,2+(i-1)*3) StressNode(j,4+(i-1)*3) 0
#                  StressNode(j,4+(i-1)*3) StressNode(j,3+(i-1)*3) 0
#                 0 0 pois*(StressNode(j,2+(i-1)*3)+StressNode(j,3+(i-1)*3))]
#             principal = eig(s)
#             Tresca(j,i) = max(principal)-min(principal)
#         end
#     end
# end
       
# StrainOut = zeros(nelem,13)
# StrainOut(1:nelem,1) = elnodes(:,1)
# StrainOut(:,2:7) = strain(:,1:6)
# StrainOut(:,8:10) = strain(:,10:12)
# StrainOut(:,11:13) = strain(:,7:9)

# SupReac = zeros(ndispl,3)
# SupReac(:,1:2) = displ(:,1:2)
# SupReac(:,3) = Pb

# fid = fopen(file_out,'w')
# fprintf(fid,'OUTPUT OF MATLAB Q4 SMALL STRAIN FEM IMPLEMENTATION \n')
# fprintf(fid,'\n')
# fprintf(fid,'          DISPLACEMENTS \n')
# fprintf(fid,'********************************* \n')
# fprintf(fid,'  Node      U1           U2 \n')
# fprintf(fid,'********************************* \n')
# fprintf(fid,'#5d #13.5e #13.5e \n',Uoutput')

# fprintf(fid,'\n')
# fprintf(fid,'                ELEMENT STRESSES \n')
# fprintf(fid,['*********************************************************', ...
#              '*********************************************************', ...
#              '*********************************************** \n'])
# fprintf(fid,['Element  S11_G1       S22_G1       S12_G1      S11_G2   ', ...
#              '    S22_G2       S12_G2       S11_G3       S22_G3       ', ...
#              ' S12_G3      S11_G4       S22_G4       S12_G4 \n'])
# fprintf(fid,['*********************************************************', ...
#              '*********************************************************', ...
#              '*********************************************** \n'])
# fprintf(fid,['#5d #12.4e #12.4e #12.4e #12.4e #12.4e #12.4e #12.4e ', ...
#              '#12.4e #12.4e #12.4e #12.4e #12.4e \n'],StressOut')

# fprintf(fid,'\n')
# fprintf(fid,'                ELEMENT STRAINS \n')
# fprintf(fid,['*********************************************************', ...
#              '*********************************************************', ...
#              '*********************************************** \n'])
# fprintf(fid,['Element  ES11_G1       E22_G1       E12_G1      E11_G2   ', ...
#              '    ES22_G2       E12_G2       E11_G3       E22_G3       ', ...
#              ' E12_G3      E11_G4       E22_G4       E12_G4 \n'])
# fprintf(fid,['*********************************************************', ...
#              '*********************************************************', ...
#              '*********************************************** \n'])
# fprintf(fid,['#5d #12.4e #12.4e #12.4e #12.4e #12.4e #12.4e #12.4e ', ...
#              '#12.4e #12.4e #12.4e #12.4e #12.4e \n'],StrainOut')

# fprintf(fid,'\n')
# fprintf(fid,'       SUPPORT REACTIONS \n')
# fprintf(fid,'***************************** \n')
# fprintf(fid,'  Node   Dof      Magnitude \n')
# fprintf(fid,'***************************** \n')
# fprintf(fid,'#5d #5d #17.5e \n',SupReac')
# fclose(fid)

# finish = toc
# disp(['Done writing output             : ',num2str(finish),' seconds.'])
    
# #Graphical output
# if ~visual
#     return
# end
# dpm=[U(1:2:2*nnodes) U(2:2:2*nnodes)]
# maxdispx = max(dpm(:,1))
# mindispx = min(dpm(:,1))
# maxdispy = max(dpm(:,2))
# mindispy = min(dpm(:,2))
# maxcoordx = max(coor(:,1))
# mincoordx = min(coor(:,1))
# maxcoordy = max(coor(:,2))
# mincoordy = min(coor(:,2))

# magfacx = (maxcoordx-mincoordx)/(maxdispx-mindispx)
# magfacy = (maxcoordy-mincoordy)/(maxdispy-mindispy)

# magfac = min([magfacx magfacy])

# dcoor = coor + magfac.*dpm

# disp('Click on an option in the plot menu window.')
# disp(['Please be patient in the case of large meshes. It may take a few seconds .......'])

# choice = 0
# manual = 0
# while choice ~= 12
#     choice = menu('Plot options','Undeformed shape','Deformed shape', ...
#                   'Overlay','S11 contour','S22 contour','S12 contour', ...
#                   'Von Mises contour','Tresca contour', ...
#                   'Manual colorbar scaling', ...
#                   'Automatic colorbar scaling','Set magnification factor','Quit')
#     if choice == 1
#         figure(1)
#         clf
#         hold on
#         for i=1:nelem
#             h=patch(coor(elnodes(i,[2 3 4 5 2]),1),coor(elnodes(i,[2 3 4 5 2]),2),'b')
#             set(h,'FaceAlpha',0.8)
#         end 
#         hold off
#         axis equal
#         title('Undeformed shape')
#     elseif choice == 2
#         figure(1)
#         clf
#         hold on
#         for i=1:nelem            
#             h=patch(dcoor(elnodes(i,[2 3 4 5 2]),1),dcoor(elnodes(i,[2 3 4 5 2]),2),'b')
#             set(h,'FaceAlpha',0.8)
#         end
#         hold off
#         axis equal
#         title(['Deformed shape, magnification factor = ',num2str(magfac)])
#     elseif choice == 3
#         figure(1)
#         clf
#         hold on
#         for i=1:nelem            
#             h=patch(coor(elnodes(i,[2 3 4 5 2]),1),coor(elnodes(i,[2 3 4 5 2]),2),'c')
#             set(h,'FaceAlpha',0.2)
#         end
#         for i=1:nelem            
#             h=patch(dcoor(elnodes(i,[2 3 4 5 2]),1),dcoor(elnodes(i,[2 3 4 5 2]),2),'b')
#             set(h,'FaceAlpha',0.2)
#         end
#         hold off
#         axis equal
#         title(['Deformed over undeformed shape, MagFac =',num2str(magfac)])
#     elseif choice == 4
#         figure(1)
#         clf
#         hold on
#         if manual==0
#             maxstress = max(max(StressNode(:,[2 5 8 11])))
#             minstress = min(min(StressNode(:,[2 5 8 11])))
#         end
#         for i=1:nelem
#            x=dcoor(elnodes(i,2:5),1)
#            y=dcoor(elnodes(i,2:5),2)
#            h=fill(x,y,StressNode(i,[2 5 8 11]))
#            set(h,'LineStyle','none')
#         end
#         hold off
#         caxis([minstress maxstress])
#         axis equal
#         colorbar('horiz')
#         title('Contour plot: \sigma_{11}')
#     elseif choice==5
#         figure(1)
#         clf
#         hold on
#         if manual == 0
#             maxstress = max(max(StressNode(:,[3 6 9 12])))
#             minstress = min(min(StressNode(:,[3 6 9 12])))
#         end
#         for i=1:nelem
#            x=dcoor(elnodes(i,2:5),1)
#            y=dcoor(elnodes(i,2:5),2)
#            h=fill(x,y,StressNode(i,[3 6 9 12]))
#            set(h,'LineStyle','none')
#         end
#         hold off
#         caxis([minstress maxstress])
#         axis equal
#         colorbar('horiz')
#         title('Contour plot: \sigma_{22}')
#     elseif choice==6
#         figure(1)
#         clf
#         hold on
#         if manual == 0
#             maxstress = max(max(StressNode(:,[4 7 10 13])))
#             minstress = min(min(StressNode(:,[4 7 10 13])))
#         end
#         for i=1:nelem
#            x=dcoor(elnodes(i,2:5),1)
#            y=dcoor(elnodes(i,2:5),2)
#            h=fill(x,y,StressNode(i,[4 7 10 13]))
#            set(h,'LineStyle','none')
#         end
#         hold off
#         caxis([minstress maxstress])
#         axis equal
#         colorbar('horiz')
#         title('Contour plot: \sigma_{12}')
#     elseif choice==7
#         figure(1)
#         clf
#         hold on
#         if manual == 0
#             maxstress = max(max(VonMises(:,:)))
#             minstress = min(min(VonMises(:,:)))
#         end
#         for i=1:nelem
#            x=dcoor(elnodes(i,2:5),1)
#            y=dcoor(elnodes(i,2:5),2)
#            h=fill(x,y,VonMises(i,:))
#            set(h,'LineStyle','none')
#         end
#         hold off
#         caxis([minstress maxstress])
#         axis equal
#         colorbar('horiz')
#         title('Contour plot: Von Mises')       
#       elseif choice==8
#         figure(1)
#         clf
#         hold on
#         if manual == 0
#             maxstress = max(max(Tresca(:,:)))
#             minstress = min(min(Tresca(:,:)))
#         end
#         for i=1:nelem
#            x=dcoor(elnodes(i,2:5),1)
#            y=dcoor(elnodes(i,2:5),2)
#            h=fill(x,y,Tresca(i,:))
#            set(h,'LineStyle','none')
#         end
#         hold off
#         caxis([minstress maxstress])
#         axis equal
#         colorbar('horiz')
#         title('Contour plot: Tresca')       
#     elseif choice == 9
#         oldmin = minstress
#         oldmax = maxstress
#         manual = 1
#         minstress = input('Minimum contour level ? ')
#         maxstress = input('Maximum contour level ? ')
#         caxis([minstress maxstress])
#         colorbar('horiz')
#         figure(1)
#         refresh
#     elseif choice == 10
#         manual = 0        
#     elseif choice == 11
#         magfac = input('Enter new displacement magnification factor. ')
#         dcoor = coor + magfac.*dpm
#     end
# end
        
