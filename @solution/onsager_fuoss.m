function [mobility]=onsager_fuoss(obj)
			
	% Initialize the empty variables that you will need.
	% Omega = mobility / F / z
	omega=[];
	% z_list is an unsorted list of all charges of all ions
	z_list=[];
	% conc_list is an unsorted list of all concentrations of all charge states. 
	conc_list=[];
			
	%All three share the same order
			
	%populate them.
	for i=1:length(obj.ions)
		omega=cat(1, omega, obj.ions{i}.fi_mobility./obj.F./obj.ions{i}.z);
		z_list=cat(1, z_list, obj.ions{i}.z);
		conc_list=cat(1, conc_list, (obj.concentrations(i).*obj.ions{i}.ionization_fraction(obj.pH, obj.I))');	
	end
			
	% Here is the mobility, which for some reason =mobility/F. 
	mob=omega.*abs(z_list);
			
	n_states=length(omega);
			
	% potential is the (chemical?) potential of each ion.
	potential=conc_list.*z_list.^2./obj.I/2; %total potential
	
	% initialize and populate the h matrix. 
    h=zeros(n_states);
	for j=1:n_states
	    for i=1:n_states
	   	h(j,i)=potential(i)*omega(i)/(omega(i)+omega(j));
	    end
	end
	
	d=diag(sum(h,2));
	h=h+d;
	II=eye(n_states);
			 
    B=2*h-II;
			
    r(:,1) = (z_list-sum(z_list.*potential) / sum(potential.*abs(z_list)./mob) * (abs(z_list)./mob))'; %check for absolute signs 

    for i=2:6                                        
    	r(:,i)=B*r(:,i-1);
    end
			
	c=[0.2929 -0.3536 0.0884 -0.0442 0.0276 -0.0193];
	%coefficients in onsager-fuoss paper

    factor=(c*r')';
    mob_new=omega-(0.78420*z_list.*factor.*omega+31.410e-9/obj.F).*sqrt(obj.I)./(1+1.5*sqrt(obj.I));
	mob_new=(mob_new.*z_list*obj.F)';
		
		
		
	index=1;
	mobility=cell(1, length(obj.ions));
	for i=1:length(obj.ions)
		mobility{i}=mob_new(index:(index+length(obj.ions{i}.z)-1));
		index=index+length(obj.ions{i}.z);
	end
end