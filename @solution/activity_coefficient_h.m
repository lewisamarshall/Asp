function gam_h=activity_coefficient_h(obj, I)
	if ~exist('I', 'var')
		I=obj.I;
	end
			
	% There are two coefficients that are used repeatedly.
	% Specified in Bahga. 
	A=obj.Adh*sqrt(I)/(1+obj.aD*sqrt(I));
	B=obj.Adh*0.1*I; %altered to match spresso code, Adh may not belong
			
	% Use them to calculate the activity coefficients. 
	% These coefficients are for z=+-1, for H+ and OH-
	gam_h=B-A;
	gam_h=10^gam_h;
end