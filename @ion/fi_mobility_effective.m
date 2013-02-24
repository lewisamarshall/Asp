function fi_mobility_effective=get_fi_mobility_effective(obj, I)
	% If only an ionic strength is specified, use the Robinson-Stokes 
	% correction to calculate a new fully ionized mobility. 
			
	% If a solution object is supplied, use the full onsager fouss correction. 
			
	if isnumeric(I) && I>=0
		% Robinson-Stokes correction. Currently using the ionic strength where Bahga 2010
		% uses twice the ionic strength. This appears to work, and follows the SPRESSO implimentation. 
		% Likely typo in paper. 
		A=0.2297;
		B=31.410e-9;
		fi_mobility_effective=obj.fi_mobility-(A.*obj.fi_mobility+B.*sign(obj.z))*sqrt(I)/(1+1.5*sqrt(I));
				
	elseif isa(I, 'solution')
		TMP=I.add_ion(obj, 0);
		fi_mobility_effective=TMP.get_factor;
		fi_mobility_effective=fi_mobility_effective(end-length(obj.z)+1:end);
	else
		error('Ionic strength must be specified as either a positive scalar value or a solution object.')
	end
end