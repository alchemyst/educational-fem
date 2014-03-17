%function quad4(filename, file_out, visual)

clear all
filename = 'sample.inp';
file_out = 'sample.out';
visual = false;

% Open input file
% if ~exist('filename', 'var')
%     filename = input('Input filename ? ','s');
% end
% if ~exist('file_out', 'var')
%     file_out = input('Output filename ? ','s');
% end
% if ~exist('visual', 'var')
%     visual = false;
% end

fid = fopen(filename,'r');

tic;
% Read number of nodes
dummy = fgetl(fid);
nnodes = fscanf(fid,'%d \n',1);

%Read nodal coordinates
dummy = fgetl(fid);
ndcoor = fscanf(fid,'%f \n',[3,nnodes])';

% Store node numbers in nodes
nodes = round(ndcoor(:,1));

% Store coordinates in coord
coor = ndcoor(:,2:3);
clear ndcoor

% Read number of elements
dummy = fgetl(fid);
nelem = fscanf(fid,'%d \n',1);

% Read if elements are plane stress or plain strain
dummy = fgetl(fid);
plane = fscanf(fid,'%1d \n',1);

% Read element number and element connectivity 
dummy = fgetl(fid);
elnodes = fscanf(fid,'%d \n',[5,nelem])';

% Read material constants
dummy = fgetl(fid);
elas = fscanf(fid,'%f',1);
pois = fscanf(fid,'%f',1);

% Read element thickness
t = fscanf(fid,'%f \n',1);

% Read number of prescribed displacements
dummy = fgetl(fid);
ndispl = fscanf(fid,'%d \n',1);

% Read prescribed displacements;
dummy = fgetl(fid);
displ = fscanf(fid,'%f \n',[3,ndispl])';
displ = sortrows(displ);

% Read number of nodal loads
dummy = fgetl(fid);
ncload = fscanf(fid,'%d \n',1);

% Read nodal loads
dummy = fgetl(fid);
cload = fscanf(fid,'%f \n',[3,ncload])';
fclose(fid);
finish = toc;
disp(['Done reading input file         : ',num2str(finish),' seconds'])

tic;
% Find specified (sdof) and unknown (udof) degrees of freedom
dof = ones(nnodes*2,1);
for i=1:ndispl
    pos=(find(nodes==displ(i,1))-1)*2+displ(i,2);
    dof(pos,1)=0;
    Ub(i,1)=displ(i,3);
end
sdof=find(dof==0);
udof=find(dof~=0);

% Construct D matrix
% If plane strain
if plane==0
    e = elas/(1-pois^2);
    nu = pois/(1-pois);
else
    e = elas;
    nu = pois;
end

c = e/(1-nu^2);
E = zeros(3,3);
E(1,1) = c;
E(2,2) = c;
E(1,2) = c*nu;
E(2,1) = E(1,2);
E(3,3) = 0.5*c*(1-nu);

%Initialize global stiffness matrix;
row_vec   = zeros(64*nelem,1);
col_vec   = zeros(64*nelem,1);
stiff_vec = zeros(64*nelem,1);
 
%Initialize global load vector
p_global = zeros(2*nnodes,1);

%Initialize B_ALL and DB_ALL
B_ALL  = zeros(12*nelem,8);
EB_ALL = zeros(12*nelem,8);

Gauss_pos = 1.0/sqrt(3.0);
pos_vec = 0;
%Main loop over elements. Compute k_elem and assemble
B     = zeros(3,8);
D     = zeros(4,8);
dNdxi = zeros(2,4);
B3    = zeros(3,4);
invJ  = zeros(2,2);
pg    = zeros(8,1);
for i=1:nelem
    % Find coordinates of element nodes
    x = coor(elnodes(i,2:5),1);
    y = coor(elnodes(i,2:5),2);
    XY = [x y];
    % Numerical integration loops    
    k_elem  = zeros(8,8);
    Gausspt = 0;
    for jGauss = 1:2
        for iGauss = 1:2
            xi  = (-1)^iGauss*Gauss_pos;
            eta = (-1)^jGauss*Gauss_pos;
            dNdxi(1,1) = 0.25*(eta - 1.0);
            dNdxi(1,2) = 0.25*(1.0 - eta);
            dNdxi(1,3) = 0.25*(1.0 + eta);
            dNdxi(1,4) = 0.25*(-1.0 - eta);
            dNdxi(2,1) = 0.25*(xi - 1.0);
            dNdxi(2,2) = 0.25*(-1.0 - xi);
            dNdxi(2,3) = 0.25*(1.0 + xi);
            dNdxi(2,4) = 0.25*(1.0 - xi);
            D(1:2,1:2:7) = dNdxi;
            D(3:4,2:2:8) = dNdxi;
            J = dNdxi*XY;
            detJ = J(1,1)*J(2,2) - J(2,1)*J(1,2);
            invJ(1,1) =  J(2,2)/detJ;
            invJ(1,2) = -J(1,2)/detJ;
            invJ(2,1) = -J(2,1)/detJ;
            invJ(2,2) =  J(1,1)/detJ;
            B2(1,1) = invJ(1,1);
            B2(1,2) = invJ(1,2);
            B2(2,3) = invJ(2,1);
            B2(2,4) = invJ(2,2);
            B2(3,1) = invJ(2,1);
            B2(3,2) = invJ(2,2);
            B2(3,3) = invJ(1,1);
            B2(3,4) = invJ(1,2);
            B = B2*D;
            k_elem = k_elem + t*B'*E*B*detJ;
            B_ALL((i-1)*12+Gausspt*3+1:(i-1)*12+Gausspt*3+3,:)=B;
            EB_ALL((i-1)*12+Gausspt*3+1:(i-1)*12+Gausspt*3+3,:)=E*B;
            Gausspt = Gausspt + 1;
        end
    end
% Assemble k_elem into k_global
    p     = elnodes(i,2:5);
    pg(1) = 2*p(1)-1;
    pg(2) = 2*p(1);
    pg(3) = 2*p(2)-1;
    pg(4) = 2*p(2);
    pg(5) = 2*p(3)-1;
    pg(6) = 2*p(3);
    pg(7) = 2*p(4)-1;
    pg(8) = 2*p(4);
     for ii=1:8
         for jj=1:8
            % if (pg(jj)>=pg(ii))
                  pos_vec = pos_vec + 1;
                  row_vec(pos_vec)   = pg(ii);
                  col_vec(pos_vec)   = pg(jj);
                  stiff_vec(pos_vec) = k_elem(ii,jj);
            % end
         end
     end
end;  % End of main loop over elements
k_global = sparse(row_vec,col_vec,stiff_vec,2*nnodes,2*nnodes);
%clear pos_vec row_vec stiff_vec

% Add nodal loads to global load vector
for i=1:ncload;
    p = find(nodes==cload(i,1));
    pos = (p-1)*2+cload(i,2);
    p_global(pos,1) = p_global(pos,1)+cload(i,3);
end

%Partition k_global into Kaa, Kab, Kba, Kbb
Kaa = k_global(udof,udof);
Kab = k_global(udof,sdof);
Kbb = k_global(sdof,sdof);
%clear k_global

% Modify RHS to include specified displacements
Pa  = p_global(udof) - Kab*Ub;

finish = toc;
disp(['Done assembling stiffness matrix: ',num2str(finish),' seconds.'])

tic;
%Solve unknown dof's
Kaa = chol(Kaa);

Ua  = Kaa'\Pa;
Ua  = Kaa\Ua;
finish = toc;

disp(['Done solving system             : ',num2str(finish),' seconds.'])

tic;

% Solve support reactions
Pb = Kab'*Ua + Kbb*Ub;

% Sort Ua and Ub into A
U(udof,1) = Ua;
U(sdof,1) = Ub;

%Compute element strains and stresses
for i=1:nelem
    p  = elnodes(i,2:5);
    pg = 2*[p(1)-0.5 p(1) p(2)-0.5 p(2) p(3)-0.5 p(3) p(4)-0.5 p(4)];
    strain(i,:) = (B_ALL((i-1)*12+1:12*i,:)*U(pg))';
    stress(i,:) = (EB_ALL((i-1)*12+1:12*i,:)*U(pg))';
end

% Write output to file
Uoutput = zeros(nnodes,3);
Uoutput(:,1) = nodes;
Uoutput(:,2) = U(1:2:nnodes*2);
Uoutput(:,3) = U(2:2:nnodes*2);

StressOut = zeros(nelem,13);
StressOut(:,1)     = elnodes(:,1);
StressOut(:,2:7)   = stress(:,1:6);
StressOut(:,8:10)  = stress(:,10:12);
StressOut(:,11:13) = stress(:,7:9);

StressNode(1:nelem,1) = elnodes(:,1);
StressNode(:,2) = (1+sqrt(3)/2)*StressOut(:,2) - 0.5*(StressOut(:,5) + ...
                  StressOut(:,11)) + (1-sqrt(3)/2)*StressOut(:,8);
StressNode(:,5) = (1+sqrt(3)/2)*StressOut(:,5) - 0.5*(StressOut(:,2) + ...
                  StressOut(:,8)) + (1-sqrt(3)/2)*StressOut(:,11);
StressNode(:,8) = (1+sqrt(3)/2)*StressOut(:,8) - 0.5*(StressOut(:,5) + ...
                  StressOut(:,11)) + (1-sqrt(3)/2)*StressOut(:,2);
StressNode(:,11) = (1+sqrt(3)/2)*StressOut(:,11) - 0.5*(StressOut(:,2) + ...
                  StressOut(:,8)) + (1-sqrt(3)/2)*StressOut(:,5);

StressNode(:,3) = (1+sqrt(3)/2)*StressOut(:,3) - 0.5*(StressOut(:,6) + ...
                  StressOut(:,12)) + (1-sqrt(3)/2)*StressOut(:,9);
StressNode(:,6) = (1+sqrt(3)/2)*StressOut(:,6) - 0.5*(StressOut(:,3) + ...
                  StressOut(:,9)) + (1-sqrt(3)/2)*StressOut(:,12);
StressNode(:,9) = (1+sqrt(3)/2)*StressOut(:,9) - 0.5*(StressOut(:,6) + ...
                  StressOut(:,12)) + (1-sqrt(3)/2)*StressOut(:,3);
StressNode(:,12) = (1+sqrt(3)/2)*StressOut(:,12) - 0.5*(StressOut(:,3) + ...
                  StressOut(:,9)) + (1-sqrt(3)/2)*StressOut(:,6);

StressNode(:,4) = (1+sqrt(3)/2)*StressOut(:,4) - 0.5*(StressOut(:,7) + ...
                  StressOut(:,13)) + (1-sqrt(3)/2)*StressOut(:,10);
StressNode(:,7) = (1+sqrt(3)/2)*StressOut(:,7) - 0.5*(StressOut(:,4) + ...
                  StressOut(:,10)) + (1-sqrt(3)/2)*StressOut(:,13);
StressNode(:,10) = (1+sqrt(3)/2)*StressOut(:,10) - 0.5*(StressOut(:,7) + ...
                  StressOut(:,13)) + (1-sqrt(3)/2)*StressOut(:,4);
StressNode(:,13) = (1+sqrt(3)/2)*StressOut(:,13) - 0.5*(StressOut(:,4) + ...
                  StressOut(:,10)) + (1-sqrt(3)/2)*StressOut(:,7);

VonMises = zeros(nelem,4);
if (plane==1)
    for i=1:4
        VonMises(:,i) = StressNode(:,2+(i-1)*3).^2 - ...
                StressNode(:,2+(i-1)*3).*StressNode(:,3+(i-1)*3) + ...
                StressNode(:,3+(i-1)*3).^2 + 3*StressNode(:,4+(i-1)*3).^2;
        VonMises(:,i) = VonMises(:,i).^0.5;
    end
else
    for i=1:4
        VonMises(:,i) = (1-pois+pois^2)*(StressNode(:,2+(i-1)*3).^2 + ...
                        StressNode(:,3+(i-1)*3).^2) - (1+pois-pois^2)* ...
                        StressNode(:,2+(i-1)*3).*StressNode(:,3+(i-1)*3)+...
                        3*StressNode(:,4+(i-1)*3).^2;
        VonMises(:,i) = VonMises(:,i).^0.5;
    end
end    

Tresca = zeros(nelem,4);
if plane ==1
    for j=1:nelem
        for i=1:4
            s = [StressNode(j,2+(i-1)*3) StressNode(j,4+(i-1)*3) 0
                 StressNode(j,4+(i-1)*3) StressNode(j,3+(i-1)*3) 0
                0 0 0];
            principal = eig(s);
            Tresca(j,i) = max(principal)-min(principal);
        end
    end
else
    for j=1:nelem
        for i=1:4
            s = [StressNode(j,2+(i-1)*3) StressNode(j,4+(i-1)*3) 0
                 StressNode(j,4+(i-1)*3) StressNode(j,3+(i-1)*3) 0
                0 0 pois*(StressNode(j,2+(i-1)*3)+StressNode(j,3+(i-1)*3))];
            principal = eig(s);
            Tresca(j,i) = max(principal)-min(principal);
        end
    end
end
       
StrainOut = zeros(nelem,13);
StrainOut(1:nelem,1) = elnodes(:,1);
StrainOut(:,2:7) = strain(:,1:6);
StrainOut(:,8:10) = strain(:,10:12);
StrainOut(:,11:13) = strain(:,7:9);

SupReac = zeros(ndispl,3);
SupReac(:,1:2) = displ(:,1:2);
SupReac(:,3) = Pb;

fid = fopen(file_out,'w');
fprintf(fid,'OUTPUT OF MATLAB Q4 SMALL STRAIN FEM IMPLEMENTATION \n');
fprintf(fid,'\n');
fprintf(fid,'          DISPLACEMENTS \n');
fprintf(fid,'********************************* \n');
fprintf(fid,'  Node      U1           U2 \n');
fprintf(fid,'********************************* \n');
fprintf(fid,'%5d %13.5e %13.5e \n',Uoutput');

fprintf(fid,'\n');
fprintf(fid,'                ELEMENT STRESSES \n');
fprintf(fid,['*********************************************************', ...
             '*********************************************************', ...
             '*********************************************** \n']);
fprintf(fid,['Element  S11_G1       S22_G1       S12_G1      S11_G2   ', ...
             '    S22_G2       S12_G2       S11_G3       S22_G3       ', ...
             ' S12_G3      S11_G4       S22_G4       S12_G4 \n']);
fprintf(fid,['*********************************************************', ...
             '*********************************************************', ...
             '*********************************************** \n']);
fprintf(fid,['%5d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e ', ...
             '%12.4e %12.4e %12.4e %12.4e %12.4e \n'],StressOut');

fprintf(fid,'\n');
fprintf(fid,'                ELEMENT STRAINS \n');
fprintf(fid,['*********************************************************', ...
             '*********************************************************', ...
             '*********************************************** \n']);
fprintf(fid,['Element  ES11_G1       E22_G1       E12_G1      E11_G2   ', ...
             '    ES22_G2       E12_G2       E11_G3       E22_G3       ', ...
             ' E12_G3      E11_G4       E22_G4       E12_G4 \n']);
fprintf(fid,['*********************************************************', ...
             '*********************************************************', ...
             '*********************************************** \n']);
fprintf(fid,['%5d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e ', ...
             '%12.4e %12.4e %12.4e %12.4e %12.4e \n'],StrainOut');

fprintf(fid,'\n');
fprintf(fid,'       SUPPORT REACTIONS \n');
fprintf(fid,'***************************** \n');
fprintf(fid,'  Node   Dof      Magnitude \n');
fprintf(fid,'***************************** \n');
fprintf(fid,'%5d %5d %17.5e \n',SupReac');
fclose(fid);

finish = toc;
disp(['Done writing output             : ',num2str(finish),' seconds.'])
    
%Graphical output
if ~visual
    return
end
dpm=[U(1:2:2*nnodes) U(2:2:2*nnodes)];
maxdispx = max(dpm(:,1));
mindispx = min(dpm(:,1));
maxdispy = max(dpm(:,2));
mindispy = min(dpm(:,2));
maxcoordx = max(coor(:,1));
mincoordx = min(coor(:,1));
maxcoordy = max(coor(:,2));
mincoordy = min(coor(:,2));

magfacx = (maxcoordx-mincoordx)/(maxdispx-mindispx);
magfacy = (maxcoordy-mincoordy)/(maxdispy-mindispy);

magfac = min([magfacx magfacy]);

dcoor = coor + magfac.*dpm;

disp('Click on an option in the plot menu window.')
disp(['Please be patient in the case of large meshes. It may take a few seconds .......'])

choice = 0;
manual = 0;
while choice ~= 12
    choice = menu('Plot options','Undeformed shape','Deformed shape', ...
                  'Overlay','S11 contour','S22 contour','S12 contour', ...
                  'Von Mises contour','Tresca contour', ...
                  'Manual colorbar scaling', ...
                  'Automatic colorbar scaling','Set magnification factor','Quit');
    if choice == 1
        figure(1)
        clf
        hold on
        for i=1:nelem
            h=patch(coor(elnodes(i,[2 3 4 5 2]),1),coor(elnodes(i,[2 3 4 5 2]),2),'b');
            set(h,'FaceAlpha',0.8)
        end 
        hold off
        axis equal
        title('Undeformed shape')
    elseif choice == 2
        figure(1)
        clf
        hold on
        for i=1:nelem            
            h=patch(dcoor(elnodes(i,[2 3 4 5 2]),1),dcoor(elnodes(i,[2 3 4 5 2]),2),'b');
            set(h,'FaceAlpha',0.8)
        end
        hold off
        axis equal
        title(['Deformed shape, magnification factor = ',num2str(magfac)])
    elseif choice == 3
        figure(1)
        clf
        hold on
        for i=1:nelem            
            h=patch(coor(elnodes(i,[2 3 4 5 2]),1),coor(elnodes(i,[2 3 4 5 2]),2),'c');
            set(h,'FaceAlpha',0.2)
        end
        for i=1:nelem            
            h=patch(dcoor(elnodes(i,[2 3 4 5 2]),1),dcoor(elnodes(i,[2 3 4 5 2]),2),'b');
            set(h,'FaceAlpha',0.2);
        end
        hold off
        axis equal
        title(['Deformed over undeformed shape, MagFac =',num2str(magfac)])
    elseif choice == 4
        figure(1)
        clf
        hold on
        if manual==0
            maxstress = max(max(StressNode(:,[2 5 8 11])));
            minstress = min(min(StressNode(:,[2 5 8 11])));
        end
        for i=1:nelem
           x=dcoor(elnodes(i,2:5),1);
           y=dcoor(elnodes(i,2:5),2);
           h=fill(x,y,StressNode(i,[2 5 8 11]));
           set(h,'LineStyle','none')
        end
        hold off
        caxis([minstress maxstress]);
        axis equal
        colorbar('horiz')
        title('Contour plot: \sigma_{11}')
    elseif choice==5
        figure(1)
        clf
        hold on
        if manual == 0
            maxstress = max(max(StressNode(:,[3 6 9 12])));
            minstress = min(min(StressNode(:,[3 6 9 12])));
        end
        for i=1:nelem
           x=dcoor(elnodes(i,2:5),1);
           y=dcoor(elnodes(i,2:5),2);
           h=fill(x,y,StressNode(i,[3 6 9 12]));
           set(h,'LineStyle','none')
        end
        hold off
        caxis([minstress maxstress]);
        axis equal
        colorbar('horiz')
        title('Contour plot: \sigma_{22}')
    elseif choice==6
        figure(1)
        clf
        hold on
        if manual == 0
            maxstress = max(max(StressNode(:,[4 7 10 13])));
            minstress = min(min(StressNode(:,[4 7 10 13])));
        end
        for i=1:nelem
           x=dcoor(elnodes(i,2:5),1);
           y=dcoor(elnodes(i,2:5),2);
           h=fill(x,y,StressNode(i,[4 7 10 13]));
           set(h,'LineStyle','none')
        end
        hold off
        caxis([minstress maxstress]);
        axis equal
        colorbar('horiz')
        title('Contour plot: \sigma_{12}')
    elseif choice==7
        figure(1)
        clf
        hold on
        if manual == 0
            maxstress = max(max(VonMises(:,:)));
            minstress = min(min(VonMises(:,:)));
        end
        for i=1:nelem
           x=dcoor(elnodes(i,2:5),1);
           y=dcoor(elnodes(i,2:5),2);
           h=fill(x,y,VonMises(i,:));
           set(h,'LineStyle','none')
        end
        hold off
        caxis([minstress maxstress]);
        axis equal
        colorbar('horiz')
        title('Contour plot: Von Mises')       
      elseif choice==8
        figure(1)
        clf
        hold on
        if manual == 0
            maxstress = max(max(Tresca(:,:)));
            minstress = min(min(Tresca(:,:)));
        end
        for i=1:nelem
           x=dcoor(elnodes(i,2:5),1);
           y=dcoor(elnodes(i,2:5),2);
           h=fill(x,y,Tresca(i,:));
           set(h,'LineStyle','none')
        end
        hold off
        caxis([minstress maxstress]);
        axis equal
        colorbar('horiz')
        title('Contour plot: Tresca')       
    elseif choice == 9
        oldmin = minstress;
        oldmax = maxstress;
        manual = 1;
        minstress = input('Minimum contour level ? ');
        maxstress = input('Maximum contour level ? ');
        caxis([minstress maxstress]);
        colorbar('horiz')
        figure(1)
        refresh
    elseif choice == 10
        manual = 0;        
    elseif choice == 11
        magfac = input('Enter new displacement magnification factor. ');
        dcoor = coor + magfac.*dpm;
    end
end
        
