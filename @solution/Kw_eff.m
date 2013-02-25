function Kw_eff=Kw_eff(obj, I)
	% Gets the effective water dissociation constant
	% due to activity corrections to H+ and OH-. 
	if ~exist('I', 'var')
		I=obj.I;
	end
			
	gam_h=obj.H.activity_coefficient(I,1);
	Kw_eff=obj.Kw/gam_h^2;
end