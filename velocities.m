


mu=398600.4418; %standard gravitational parameter in km^3/s^2
Re=6371; %Earth Radius in km
altmin=200; %orbital altitude in km
altmax=400;

velmin=sqrt(mu/(altmin+Re));
velmax=sqrt(mu/(altmax+Re));

vel=linspace(velmin,velmax,25); %magnitude of the orbital velocity in m/s


for i=1:length(vel)

velocity=[vel(i); 0; 0;]; 

cyclev(velocity);

system('./ptetra>>/dev/null');
printf('%d%% complete...\r',(i/length(vel)*100))

fid=fopen('momout.txt' ,'r') ;

drag(i,1)=fscanf(fid,'Drag Force:\nx=   %f\n');
drag(i,2)=fscanf(fid,'y=   %f\n');
drag(i,3)=fscanf(fid,'z=   %f\n');

torq(i,1)=fscanf(fid,'Torques:\nx=   %f\n');
torq(i,2)=fscanf(fid,'y=   %f\n');
torq(i,3)=fscanf(fid,'z=   %f\n');
fflush(stdout);
end

for i=1:length(vel)
	d(i)=norm(drag(i,:));
end

plot(vel,d)
title('Drag Force as a Function of Velocity')
xlabel('Orbital Velocity (m/s along x axis)')
ylabel('Drag Force (Newtons)')

save('velvsdrag.mat','drag','torq','vel');

