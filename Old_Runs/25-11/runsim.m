

mu=398600.4418; %standard gravitational parameter in km^3/s^2
Re=6371; %Earth Radius in km
alt=300; %orbital altitude in km

vel=sqrt(mu/(alt+Re))*1000; %magnitude of the orbital velocity in m/s

pitch=linspace(-pi,pi,100);

printf('\n Simulation Start \n\n Velocity=%f\n',vel);

for i=1:length(pitch)

velocity=[cos(pitch(i)) 0 -sin(pitch(i)); 0 1 0; sin(pitch(i)) 0 cos(pitch(i))]*[vel; 0; 0;]; 

cycle(velocity);

system('./ptetra>>/dev/null');
printf('%d%% complete...\n',(i/length(pitch)*100))

fid=fopen('momout.txt' ,'r') ;

drag(i,1)=fscanf(fid,'Drag Force:\nx=   %f\n');
drag(i,2)=fscanf(fid,'y=   %f\n');
drag(i,3)=fscanf(fid,'z=   %f\n');

torq(i,1)=fscanf(fid,'Torques:\nx=   %f\n');
torq(i,2)=fscanf(fid,'y=   %f\n');
torq(i,3)=fscanf(fid,'z=   %f\n');
fflush(stdout);
end

for i=1:length(pitch)
	d(i)=norm(drag(i,:));
end

plot(pitch,d)

save('octaveout.mat','drag','torq','pitch');

