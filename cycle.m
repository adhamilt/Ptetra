function cycle(velocity)


fid=fopen('pictetra.dat','r');

A=char(fread(fid));
fclose(fid);

%find the location of each relevant variable
for i=1:length(A)
	if(strncmp(A(i:end),'vixyz=',6))
		vixyz=i;
		break;
	end
end

for i=1:length(A)
	if(strncmp(A(i:end),'vexyz=',6))
		vexyz=i;
		break;
	end
end



for i=vixyz:length(A)
	if(strncmp(A(i),char(10),1))
		vixyzend=i;
		break;
	end
end

for i=vexyz:length(A)
	if(strncmp(A(i),char(10),1))
		vexyzend=i;
		break;
	end
end


A=[A(1:vexyz-1); sprintf('vexyz= %.2e %.2e %.2e',velocity)'; A(vexyzend:vixyz-1); sprintf('vixyz= %.2e %.2e %.2e',velocity)'; A(vixyzend:end)];

fid=fopen('pictetra.dat','w');
fprintf(fid,'%s',A);
fclose(fid);
endfunction


