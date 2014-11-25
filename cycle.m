function cycle(velocity,O,N2,O2,He,Ar,H,N,T,Tex)

if nargin<10
Tex=1016;
end

if nargin<9
T=1E-3;
end

if nargin<8
N=1.332e7;
end

if nargin<7
H=1.555e5;
end

if nargin<6
Ar=2.686e6;
end

if nargin<5
He=1.187e7;
end

if nargin<4
O2=1.945e8;
end

if nargin<3
N2=3.121e9;
end

if nargin<2
O=3.852e9;
end

if nargin<1
velocity=[7E3 0 0];
end

V=0.016128; %volume of the simulation domain in m^3

ni=[O N2 O2 He Ar H N].*1E3*V; %convert densities to #/m^3

m=[15.9999 2*14.007 2*15.9999 4.0026 39.948 1.0079 14.007];




fid=fopen('pictetra.dat');

out=[];


line=1;

while(line==1)
	
	A=fgetl(fid);
	if(strncmp(A,' nipop=',7))
		out=[out ' nipop=7' char(10)];
	else
		out=[out A char(10)];
	end
	
	if(strncmp(A,'$begin plasmaparameters',23));
		line=0;
		break;
	end
end

out=[out sprintf('  ne=0.\n  fromne=-999.\n  te=1.0e-3\n  vexyz=  1.00e+04 0.00e+00 0.00e+00\n')];

Temp=T;

for i=1:7
	if(i>3)
		Temp=Tex;
	end
	
out=[out sprintf('  mi=%.2e\n  qi=0.\n  ni=%.2e\n  fromni=-999.\n  ti=%.2e\n  vixyz= %.2e %.2e %.2e\n', m(i), ni(i), Temp, velocity)];
end



while(line==0)
	A=fgetl(fid);
	if(strncmp(A,'$end plasmaparam',16))
	out=[out A char(10)];
	line=1;
	break;
	end
end

out=[out fread(fid)'];

fclose(fid);
fid=fopen('pictetra.dat','w+');
fprintf(fid,'%s',out);
fclose(fid);

endfunction
