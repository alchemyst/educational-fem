h      = input('Height of rectangle [mm] ? ');
l      = input('Length of rectangle [mm] ? ');
m      = input('Increments along height ? ');
n      = input('Increments along length ? ');
option = input('0:Plane strain or 1:Plane stress ? ');
E      = input('Elasticity modulus [MPa] ? ');
nu     = input('Poisson''s ratio ? ');
t      = input('Element thickness ? ');
V      = input(['Total shear force at beam tip (positive for upward ' ...
                'load) [N] ? ']);
filen  = input('Write input to which file ? ','s');

node = 0;
deltx = l/n;
delty = h/m;
for i=1:n+1;
    for j=1:m+1;
        node = node + 1;
        x = (i-1)*deltx;
        y = (j-1)*delty - h/2;
        coord(node,:)=[node x y];
    end
end

el = 0;
for i=1:n
    for j=1:m
        elnode1 = (1+m)*(i-1) + j;
        elnode2 = (1+m)*(i-1) +j + 1;
        elnode3 = (1+m)*i + j;
        elnode4 = (1+m)*i + j + 1;
        el1 = el+1;
        elnode(el1,:) = [el1 elnode1 elnode3 elnode4 elnode2];
        el=el+1;    
    end
end

%Compute nodal forces due to parabolic distributed shear force at beam end
delty = h / m;
forces = zeros(m+1,1);
for i=1:m;
    y1 = -h/2 + (i-1)*delty;
    y2 = -h/2 + i*delty;
    F1 = V/delty*(0.75/h*(y2^2-y1^2)-1.5/h^3*(y2^4-y1^4)-1.5/h*y1*(y2-y1) + ...
         2/h^3*y1*(y2^3-y1^3));
    F2 = V/delty*(-0.75/h*(y2^2-y1^2)+1.5/h^3*(y2^4-y1^4)+1.5/h*y2*(y2-y1) - ...
         2/h^3*y2*(y2^3-y1^3)); 
    forces(i:i+1) = forces(i:i+1) + [F1 F2]'; 
 end

fid=fopen(filen,'w')
fprintf(fid,'Number_of_nodes \n');
fprintf(fid,'%7d \n',node);
fprintf(fid,'Nodal_coordinates \n');
fprintf(fid,'%7d %20.15e %20.15e \n',coord');
fprintf(fid,'Number_of_elements \n');
fprintf(fid,'%7d \n',el);
fprintf(fid,'Plane_stress_or_strain \n');
fprintf(fid,'%7d \n',option);
fprintf(fid,'Element_connectivity \n');
fprintf(fid,'%7d %7d %7d %7d %7d \n',elnode');

fprintf(fid,'Material_properties \n');
fprintf(fid,'%20.15e %20.15e %20.15e \n',[E nu t]);
fprintf(fid,'Number_of_prescribed_displacements \n');
fprintf(fid,'%5d \n',2*m+2);
fprintf(fid,'Prescribed_displacements \n');
fprintf(fid,'%7d 1 0.0 \n',[1:m+1]);
fprintf(fid,'%7d 2 0.0 \n',[1:m+1]);
fprintf(fid,'Number_of_nodal_loads \n');
fprintf(fid,'%5d \n',m+1);
fprintf(fid,'Nodal_loads \n');
for i=1:m+1
    c_node = (m+1)*n+i;
    fprintf(fid,'%7d 2 %20.15f \n',[c_node forces(i)]);
end
fclose(fid);

figure(1)
plot(coord(:,2),coord(:,3),'x')
for i=1:node
    text(coord(i,2),coord(i,3),num2str(i));
end
axis equal