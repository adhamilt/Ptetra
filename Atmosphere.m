
year=2014; %Present Year
month=04; %April
day=01; 

alt=300;   %Altitude in km
lat=0;
long=0;

mu=398600.4418; %standard gravitational parameter in km^3/s^2
Re=6371; %Earth Radius in km

vel=sqrt(mu/(alt+Re))*1000; %magnitude of the orbital velocity in m/s

cmd= ['wget --post-data "model=msis&year=' num2str(year) '&month=' num2str(month) '&day=' num2str(day) '&time_flag=0&hour=1.5&geo_flag=0.&latitude=55.&longitude=45.&height=100.&profile=1&start=150.&stop=400.&step=1.&f10_7=&f10_7_3=&ap=&format=0&vars=05&vars=08&vars=09&vars=10&vars=11&vars=12&vars=13&vars=14&vars=15&vars=16&vars=17&linestyle=solid&charsize=&symbol=2&symsize=&yscale=Linear&xscale=Linear&imagex=640&imagey=480" http://omniweb.gsfc.nasa.gov/cgi/vitmo/vitmo_model.cgi -O MSIS.txt'];

printf('\n Obtaining MSIS Model\n');
system([cmd]);
printf('\n Done\n');

fid =fopen('MSIS.txt');

for i=1:27
	fgetl(fid);
end

for i=1:250
	s=fgetl(fid);
	[A]=sscanf(s,'  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f');
	h(i)=A(1);
	O(i)=A(2);
	N2(i)=A(3);
	O2(i)=A(4);
	rho(i)=A(5);
	T(i)=A(6);
	Tex(i)=A(7);
	He(i)=A(8);
	Ar(i)=A(9);
	H(i)=A(10);
	N(i)=A(11);
end

[val,I]=min(abs(h-alt));


pitch=linspace(-pi,pi,100);

printf('Simulating Satellite Drag:\n')

for i=1:length(pitch) 
velocity=[cos(pitch(i)) 0 -sin(pitch(i)); 0 1 0; sin(pitch(i)) 0 cos(pitch(i))]*[vel; 0; 0;]; 

cycle(velocity, O(I),N2(I),O2(I),He(I),Ar(I),H(I),N(I),T(I),Tex(I));

system('./ptetra>>/dev/null');

printf('%d%% complete...\r',(i/length(pitch)*100))

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

fclose(fid);


