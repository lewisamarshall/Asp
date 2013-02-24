classdef ion < handle
	% A class for ions dissolved in solution. 
    % Draws significantly on Bagha Electrophoresis 2010
	% "Ionic strength effects on electrophoretic focusing and separations"
	
	% The ion is defined by a name, and a set of charge states (z).
	% Each charge state must have an associated acidity constant (pKa).
	% Each charge state must have an associated fully ionized mobility (fi_mobility).
	
    properties
        name;
        z=0;
        pKa;
        fi_mobility; %Expected in m^2/V/s.  
    end
    
	% Properties only accessible by the funciton itself. 
	% These are constants and should not change.
	% Eventually, T may be  removed from the constants list.
	
    properties(Constant = true, GetAccess = 'private')
        F=9.6485E4;       	% Faraday's const.[C/mol]
        Lpm3=1000;			% Conversion from liters to m^3
		T=298;				% Temperature, in Kalvin
		% The following are constants in eqtn 6 of Bahga 2010. 
		Adh=0.5102; 		% L^1/2 / mol^1/2, approximate for room temperature
		aD=1.5*sqrt(2); 	% mol^-1/2 mol^-3/2, approximation
    end
    
    methods
        function obj=ion(name, z, pKa, fi_mobility)
			%Class constructor. 
            %Initialize the object, with checks on the form of the variables.
            if(nargin ==4 ) %
				%check that the name is a string
                if ischar(name)
                    obj.name=name;
                else
                    error('Ion name must be a string.')
                end
				
                % Check that z is a vector of integers
                if z==int8(z) & isvector(z)
                    obj.z=double(z);
                else
                    error('Charge states must be a vector of integers.')
                end
                
				% Check that the pKa is a vector of numbers of the same length as z. 
                if isvector(pKa) && length(obj.z)==length(pKa) &&isnumeric(fi_mobility)
                    obj.pKa=double(pKa);
                else
                    error ('pKas must be a numeric vector the same size as charge vector.')
                end
                
				%Check that the fully ionized mobility is a vector of numbers the same size as z. 
                if isvector(fi_mobility) && length(obj.z)==length(fi_mobility) && isnumeric(fi_mobility)
                    obj.fi_mobility=double(fi_mobility);
                else
                    error ('Fully ionized mobilities must be a numeric vector the same size as charge vector.')
                end
                
				% Force the sign of the fully ionized mobilities to match the sign of the charge. 
				% This command provides a warning, which you can suppress, with, for example, 
				% warning('off','all');
                if ~all(sign(obj.z)==sign(obj.fi_mobility))
                    obj.fi_mobility=abs(obj.fi_mobility).*double(sign(obj.z));
                    warning('Forcing fully ionized mobility signs to match charge signs.')
                end
                
            else
				% If there are not four inputs, the initialization fails.
                error ('Not enough inputs to describe an ion.')
            end
            
            % After storing the ion properties, ensure that the properties are sorted in order of charge. 
			% All other ion methods assume that the states will be sorted by charge. 
            obj=obj.z_sort(); 
        end
		
		function L=L(obj, I)
			% L products of acidity constants, for use in the pH calculating routine.
			% It can use ionic strength correction if an ionic strength is specified. Otherwise. It uses 
			% uncorrected acidity coefficients. 
			
			L=obj.z0;
			
			if ~exist('I', 'var')
				I=0;
			end
			
			Ka=obj.Ka_eff(I);
			
			index_0=find(L==0);
			L(index_0)=1;
			
			if index_0~=1
				for i=(index_0-1):1
					L(i)= L(i+1)*Ka(i);
				end
			end
			
			if index_0~=length(L)
				for i=(index_0+1):length(L)
					L(i)= L(i-1)/Ka(i-1);
				end
			end
			
		end
        
        function i_frac=ionization_fraction(obj, pH, I ,index)
        	% This function takes the ion and the pH and the ionic strength. 
			% It computes the fraction of ion in each of the ionization z states
			% in z. If index is specified, it will only return the ionization 
			% in ionization state z(index). 
			
			% If ionic strength is not specified, set it to zero. 
			if ~exist('I', 'var')
				I=0;
			end
			
			% Sanitize the pH input. 
			if ~isnumeric(pH)
				error('pH should be a number.')
			end
			
			if length(pH)~=1
				pH=pH(1); 
				warning('Ionization fraction only takes a single pH. Using pH(1).')
			end
			
			% Get the vector of products of acidity constants.
            L=obj.get_L(I);
			% Compute the concentration of H+ from the pH.
            cH=10.^(-pH)';
            
			% Calculate the denominator of the function for ionization fraction.
            i_frac_denom=1+sum(L.*bsxfun(@power, cH, obj.z),2);
            
			%Calculate the vector of ionization fractions
            i_frac=L.*bsxfun(@power, cH, obj.z)./i_frac_denom;
			
			% If index is specified, return only the ionization fraction of z(i).
			if exist('index', 'var')
				try
					i_frac=i_frac(index);
				catch
					% Send a warning if the index is meaningless. 
					% Still return the vector. 
					warning('Specified index is out of bounds.')
				end
			end
        end
        
        function m_cond=molar_conductivity(obj, pH, I)
            % This function takes the ion and the pH and computes the molar
            % conductivity. This function can take an ionic strength. 
			% Provides conducitivity in Siemens per meter per mole.
			
			%%%%%%%%%%%%%%%%
			% Does not yet use the ionic strength corrected mobilities.
			%%%%%%%%%%%%%%%%
            
			if ~exist('I', 'var')
				I=0;
			end
			
            i_frac=ionization_fraction(obj, pH, I);
            
            m_cond=sum(obj.F*obj.z.*i_frac.*obj.fi_mobility*obj.Lpm3);
            
        end
        
        function eff_mobility=effective_mobility(obj, pH, I)
            % This function takes the ion and a pH and an ionic strength. It calls the ionization
            % fraction function and uses this information to compute the
            % effective mobility of the ion. 
			
			%%%%%%%%%
			% This function uses ionic strength for ionization fraction 
			% calculation, but it does not use it for fi_mobility.
			%%%%%%%%%
			
			if ~exist('I', 'var')
				I=0;
			end
			
            for i = 1:length(pH)
                i_frac=ionization_fraction(obj, pH(i), I);
                eff_mobility(i)=sum(i_frac.*obj.fi_mobility);
            end
			
        end
        
        function obj=z_sort(obj)
			% This function sorts the charge state list (and the associated pKas)
			% by the value of the charge state (from negative to positive). This
			% is required for other methods to work correclty.
			
			% Sort the charges. Store the index order.
          	[obj.z,Index]=sort(obj.z);
			% Sort the pKas by the stored index order.
         	obj.pKa=obj.pKa(Index);
			% Sort the mobilities by the stored index order.
            obj.fi_mobility=obj.fi_mobility(Index);
			
			% This section will check each charge state to see if it is complete.
			% That is, if there is a charge state -2, there must be a charge state
			% -1. if there is a charge state +3, there must be a +2 and a +1.
			warn=0;
			for i=1:length(obj.z)
				if obj.z(i) < -1 && obj.z(i+1)~=obj.z(i)+1;
					warn=1;
				elseif obj.z(i)>1 && obj.z(i-1)~=obj.z(i)-1;
					warn=1;
				end
			end
			%Send a single warning if any charge state is missing
			if warn==1;
				warning('Charge states may be missing...')
			end		
        end

		function Ka=Ka(obj)
			% Calculates the Kas from the pKas. 
			Ka=10.^(-obj.pKa);
		end
		
		function z0=z0(obj)
			% Calls the list of charge states, but inserts 0 in the appropriate position. 
			
			z0=[0, obj.z];
			z0=sort(z0);
		end

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
				
			elseif strcmp(class(I), 'solution')
				TMP=I.add_ion(obj, 0);
				fi_mobility_effective=TMP.get_factor;
				fi_mobility_effective=fi_mobility_effective(end-length(obj.z)+1:end);
			else
				error('Ionic strength must be specified as either a positive scalar value or a solution object.')
			end
		end
		
    end %End of Methods
    
end %End of Object