function Ka_eff=Ka_eff(obj, I)
	% Uses the ionic strength correction function from
	% Dubye-Huckle theory to calculate the activity coefficients, 
	% and from this, compute the effective Ka values for the ion. 
			
	% If the ionic strength is zero, simply return the Ka's. 
	if I==0;
		Ka_eff=obj.Ka; 
		return
	end
			
	% Store the Ka values, rather than continually recalculating them.
	Ka=obj.Ka;
	% Make the effective Ka vector the same size as the Ka vector.
	Ka_eff=Ka;
			
	% Call the list of charge states, including 0. 
	% This is required b/c you need the activity of the uncharged species.
    z_list=obj.z0;
			 
	% There are two coefficients that are used repeatedly.
	% Specified in Bahga. 
	A=obj.Adh*sqrt(I)/(1+obj.aD*sqrt(I));
	B=0.1*I;
			
	% Calculate (the log of) the activity coefficients. 
	% First, the H+,
	gam_h=B-A;
	% then the ion states.
	gam_i=z_list.^2*(B-A);
			
	% Convert both to actual activity coefficients.
	gam_h=10^gam_h;
	gam_i=10.^gam_i;
			
	% For each acidity coefficient, get the effective 
	% coefficienty by multiplying by activities.
	for i=1:length(Ka_eff)
		Ka_eff(i)=Ka(i)*gam_i(i+1)/gam_i(i)/gam_h;
	end

end