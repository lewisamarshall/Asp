classdef ion
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
		fi_mobility_effective; %Intended to be filled by the solution class.
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
                if all(z==int8(z)) && isvector(z)
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
		
    end %End of Methods
    
end %End of Object