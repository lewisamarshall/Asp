function Kw_eff=Kw_eff(obj, I)
	% Gets the effective water dissociation constant
	% due to activity corrections to H+ and OH-. 
	if ~exist('I', 'var')
		I=obj.I;
	end
			
	% There are two coefficients that are used repeatedly.
	% Specified in Bahga. 
	A=obj.Adh*sqrt(I)/(1+obj.aD*sqrt(I));
	B=0.1*I;
			
	% Use them to calculate the activity coefficients. 
	% These coefficients are for z=+-1, for H+ and OH-
	gam_h=B-A;
	gam_h=10^gam_h;
	Kw_eff=obj.Kw/gam_h^2;
end