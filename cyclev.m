function cyclev(velocity)

if nargin<1
velocity=[7E3 0 0];
end


fid=fopen('pictetra.dat');

out=[];


line=1;


while(line==1)
	
	A=fgetl(fid);
	if(strncmp(A,' nipop=',7))
		out=[out ' nipop=1' char(10)];
	else
		out=[out A char(10)];
	end
	
	if(strncmp(A,'$begin plasmaparameters',23));
		line=0;
		break;
	end
end

out=[out sprintf('  ne=0.\n  fromne=-999.\n  te=1.0e-3\n  vexyz=  1.00e+04 0.00e+00 0.00e+00\n   mi=1.\n  qi=0.\n  ni=5.53e7\n  fromni=-999.\n  ti=1.0e-3\n  vixyz= %.2e %.2e %.2e\n',velocity)];

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
