classdef solution
    % SOLUTION	Create a  SOLUTION object.
	%
	% 	SOLUTION(IONS, CONCENTRATIONS) is an object representing aqueous solution
	% 		containing the ions in IONS at the concentrations specified in CONCENTRATIONS. 
	%	
	%	When a new solution is initialized, it will immediately calculate the equilibrium
	%		state, including the pH and the ionic strength (I) of the solution. These
	%		values willl be stored as permenant attributes of the object. Other solution 
	%		properties can be calculated by invoking the appropriate method. 
	%	
	% 	See also ion.
    
    properties(Constant = true, GetAccess = 'private')
        F=9.65E4;           % Faraday's const.				[C/mol]
        Rmu=8.31;           % Universal gas const. 			[J/mol*K]
        Temp=298;           % Temperature 					[K]
        Kw=1E-14;           % Water equilibrium constant	[mol^2]
        muH=362E-9;   		% Mobility of Hydronium 		[m^2/s*V]
        muOH=205E-9; 		% Mobility of Hydroxide 		[m^2/s*V]
        Lpm3=1000;          % Liters per meter^3			[]
        visc=1E-3;          % Dynamic viscosity (water) 	[Pa s]
		Adh=0.512; 			% L^1/2 / mol^1/2, approximate for room temperature
		aD=1.5*sqrt(2); 	% mol^-1/2 mol^-3/2, approximation
    end
    
    properties
        ions;				% Must be a cell of ion objects from the Asp class. 
        concentrations=0; 	% A vector of concentrations in molar.
        pH=7;				% Normal pH units. 
		I=0; 				% Expected in molar. 
    end
    
    methods
        
        function obj=solution(ions, concentrations)
			%Class Constructor

            if(nargin == 2)
				% If the object is not a cell, try to force it to be a cell array. 
				% This may be the case when there is only one ion. 
				if ~iscell(ions)
					ions=num2cell(ions);
				end
				
				% Check that all of the objects in IONS are in the ion class. 
                if isvector(ions) && all(strcmp(cellfun(@class, ions, 'UniformOutput', false), 'ion'))
                    obj.ions=ions;
                else
                    error('You must input a cell vector of ion objects. Use "help ion" for more information.')
                end
				
				% Check that CONCENTRATIONS is a numeric vector of the same length as IONS.
				% Also check that all of the concentrations are positive. 
                if isvector(concentrations) && length(ions)==length(concentrations) && isnumeric(concentrations) && all(concentrations>=0)
					% If the concentration is put in as a cell, change it to a vector. 
					if iscell(concentrations)
						concentrations=cell2mat(concentrations)
					end
					obj.concentrations=concentrations;
                else
                    error('The concentrations vector must be the same size as the ions vector.')
                end
            else % If the solution isn't specified with two arguments, throw an error. 
                error('Solutions must have a cell of ions and a cell or vector of concentrations.')
            end
            
			try
				[obj.pH, obj.I]=obj.find_equilibrium;
			catch
				warning('Could not find equilibrium with ionic strength corrections. Returning uncorrected pH and I. ')
	            obj.pH=obj.calc_pH;
				obj.I=obj.calc_I(obj.pH);
			end

        end

		function new_solution=add_ion(obj, new_ions, new_concentrations)
			% ADD_ION initializes a new solution with more ions.
			%	NEW_SOLUTION will contain all of the ions in the current solution
			%	plus a new set of ions from NEW_IONS at a new set of concentrations 
			%	from NEW_CONCENTRATIONS.
			new_solution=solution(cat(2, obj.ions, {new_ions}), cat(2, obj.concentrations, new_concentrations));
		end
		
        function I=calc_I(obj, pH, I_guess)
            % CALC_I calculates the ionic strength of a solution. 
			% 	CALC_I should only be used during solution initialization, 
			% 	after that, the equilibrium ionic strength is stored in obj.I. 
			%	pH must be supplied. If a guess for ionic strength is not supplied,
			%	ionic strength corrections will not be used. 
			
			if ~exist('I_guess', 'var')
				I_guess=0;
			end

			% Set the ionic strength to zero to start with. It will be counted for each ion.
			I=0; 
			
			% For each ion, add the contribution to ionic strength to the sum.
        	for i=1:length(obj.ions);
        		I=I+obj.concentrations(i)*sum(((obj.ions{i}.z).^2).*obj.ions{i}.ionization_fraction(pH, I_guess));
        	end
			
			% Add the ionic strength due to water dissociation. 
			I=I+obj.cH(pH) + obj.cOH(pH, I_guess);
			% Divide by 2 to get the correct value. 
        	I=I/2;
        end
        
        function pH=calc_pH(obj, I)
			% CALC_pH Finds the pH of the object. 
			%	If an ionic strength is specified, uses the corrected acidity constants. 
			%
			%	This function should be used only when finding the equilibrium state. 
			%	After that, the value should be pulled from obj.pH. 
			%
			% 	If ionic strength does not exist, assume it is zero.
			% 	This function is used to find the equilibrium state, 
			% 	so it cannot pull the ionic strength from the object.
			if ~exist('I', 'var')
				I=0; 
			end
			
			% Find the order of the polynomial. This is the maximum
			% size of the list of charge states in an ion. 
            MaxCol=-inf;
            for i=1:length(obj.concentrations)
                MaxCol=max(MaxCol, max(obj.ions{i}.z)-min(obj.ions{i}.z)+1);
            end
            
			% Set up the matrix of Ls, the multiplication
			% of acidity coefficients for each ion.
            LMat=zeros(length(obj.ions),MaxCol);
            
            for i=1:length(obj.ions)
                LMat(i,1:length(obj.ions{i}.z)+1)=obj.ions{i}.L(I);
            end
                    
            % Construct Q vector.
			Q=1;
			% Convolve each line of the L matrix. 
            for j=1:size(LMat,1)
                Q=conv(Q,LMat(j,:));
            end %for j
			%Convolve with water dissociation.
			Q=conv(Q, [-obj.Kw_eff(I) 0 1]);

            %Construct P matrix
            for i=1:length(obj.concentrations)
                z_list=obj.ions{i}.z0;
                
                tmp=zeros(1,size(LMat,2));
                tmp(1:length(z_list))=z_list;
                Mmod=LMat;     Mmod(i,:)=Mmod(i,:).*tmp;

                Pi=1;
                for kl=1:size(Mmod,1)
                    Pi=conv(Pi,Mmod(kl,:));
                end %for j
                
                Pi=conv([0 1],Pi);  % Convolve with P2
                PMat(i,:)=Pi;
            end %for i
            
            
            %Pad whichever is smaller, P or Q            
            SizeDiff=size(Q,2)-size(PMat,2);
            if SizeDiff>0
                PMat=[PMat,repmat(PMat(:,1)*0,1,SizeDiff)];
            elseif SizeDiff<0
                Q=[Q,repmat(0,1,SizeDiff)];
            end
                 
        	% Construct Polynomial
			% rewrite using repmat
            cTotRep=(obj.concentrations)'*ones(1,size(PMat,2));

            P=sum(cTotRep.*PMat,1);
            
            Poly=zeros(1,max(length(P), length(Q)));
            Poly(1:length(P))=Poly(1:length(P))+P;
            Poly(1:length(Q))=Poly(1:length(Q))+Q; %from QMat
            
            Poly=fliplr(Poly);
            
            %%%
            % Solve Polynomial for concentration (should work)
            %%%
            roo=roots(Poly);
            cH=roo(imag(roo)==0 & roo>0);
            pH=-log10(cH);
        end
		
		function Residual=equil_offset(obj,I_i)
			% EQUIL_OFFSET finds the error in the ionic strength.
			%	Takes an ionic strength, then uses it to calculate
			%	a new ionic strenght using the new equilibrum coefficents. 
			%	FIND_EQUILIBRIUM finds the root of this function.
			pH=obj.calc_pH(I_i);
			I_f=obj.calc_I(pH, I_i);
			
			Residual=(I_f-I_i);
		end
		
		function [pH, I]=find_equilibrium(obj)
			% FIND_EQUILIBRIUM finds the equilibrium ionic strength and pH. 
			%	It uses the fzero root finder to find the equilibrium
			% 	pH and ionic strength of a solution, using the ionic-
			% 	strength-adjusted activity coefficients. This function
			%	is called when the object is initialized.
			
			% Generate an initial ionic strength guess without using activity corrections
			I=obj.calc_I(obj.calc_pH);
			% Iterate to find the true ionic strength.
			
			OPTION=optimset('TolX', 1e-4);
			[I,fval,exitflag]=fzero(@obj.equil_offset, I, OPTION);

			if exitflag~=1
				error('Could not find equilibrium.')
			end
			% Use this final ionic strength to find the correct pH. 
			pH=obj.calc_pH(I);
		end
		
		function H_conductivity=H_conductivity(obj)
			% Calculates teh conductivity of H+.
			% Does not correct the mobility of the ion.
			H_conductivity=obj.F*obj.muH*obj.cH*obj.Lpm3;
		end
		
		function OH_conductivity=OH_conductivity(obj)
			% Calculates teh conductivity of OH+.
			% Does not correct the mobility of the ion.
			OH_conductivity=obj.F*obj.muOH*obj.cOH*obj.Lpm3;
		end
		
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
		
		function cH=cH(obj, pH)
			% Supplies the concentration of H+ in the solution.
			if  ~exist('pH', 'var')
				pH=obj.pH;
			end
			
			cH=10^(-pH);
		end
		
		function cOH=cOH(obj, pH, I)
			% Supplies the concentration of OH- in the solution. 
			if  ~exist('pH', 'var')
				pH=obj.pH;
			end
			if  ~exist('I', 'var')
				I=obj.I;
			end
			
			cOH=obj.Kw_eff(I)/obj.cH(pH);
		end
		
        function conductivity=conductivity(obj)
            % CONDUCTIVITY calcualtes the electrical conductivity of the solution. 
			%	Relies on molar conductivity calculations from ion and total conductivity 
			%	of H+ and OH-. 
			
            conductivity=0;
            for i=1:length(obj.concentrations)
                conductivity=conductivity+...
					obj.ions{i}.molar_conductivity(obj.pH, obj.I) * obj.concentrations(i);
            end
			
            conductivity=conductivity + obj.OH_conductivity; 	% Add contribution for hydroxyl conductivity
            conductivity=conductivity + obj.H_conductivity;		% Add contribution for hydronium conductivity
        end
        
        function T=transference(obj)
            %TRANSFERENCE Gets the fraction of charge carried by each of the ions.
            %	Should not precisely add to 1, because some charge is carried
            %	by protons and hydroxyls. 
            T=zeros(1,length(obj.ions));
            for i=1:length(T)
                T(i)=obj.ions{i}.molar_conductivity(obj.pH, obj.I).*obj.concentrations(i);
            end
            T=T/obj.conductivity;
        end
		        
        function Qi=zone_transfer(obj, vol)
            %Calculates the zone transfer charge of the solution at a given
            %volume. The volume of the solution is specified in liters. 
            Qi=obj.concentrations./obj.get_transference.*vol.*obj.F;
        end

		function Cb=buffering_capacity(obj)
			% BUFFERING_CAPACITY finds the buffering capacity of a solution. 
			%	This function generates an exact solution to the buffering
			%	capacity by finding the derivative of the pH with respect to
			%	the addition of an acid.
			
			% Find the smallest concentration in the solution. 
			c=obj.concentrations(obj.concentrations>0);
			c=0.1*min(c);
			% Add an acid insult at 1% the lowest concentration in the solution.
			new_sol=obj.add_ion(ion('Acid Insult', -1, -2, -1), c);
			% Find the slope of the pH. 
			Cb=abs(c/(obj.pH-new_sol.pH));
		end
		
		function [mob_new, factor]=get_factor(obj)
			
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
				for j=1:length(obj.ions{i}.z)
					z_list=cat(1, z_list, obj.ions{i}.z(j));
					omega=cat(1, omega, obj.ions{i}.fi_mobility(j)/obj.F/z_list(end));
				end
				conc_list=cat(1, conc_list, (obj.concentrations(i).*obj.ions{i}.ionization_fraction(obj.pH))');	
			end
			
			% Here is the mobility, which for some reason =mobility/F. 
			mob=omega.*abs(z_list);
			
			n_states=length(omega);
			mob_new=zeros(1,n_states);

    
			% mu is the (chemical?) potential of each ion.
			mu=conc_list.*z_list.^2./obj.I/2; %total potential
			
			
			for j=1:n_states
			    for i=1:n_states
			   	h(j,i)=mu(i)*omega(i)/(omega(i)+omega(j));
			    end
			end
			d=sum(h,2);
			d=diag(d);
			h=h+d;
 			II=eye(n_states,n_states);
			 
		    B=2*h-II;
			
		    r(:,1)=(z_list-sum(z_list.*mu)/sum(mu.*abs(z_list)./mob)*(abs(z_list)./mob))'; %check for absolute signs 

		    for i=2:6                                        
		    r(:,i)=B*r(:,i-1);
			end
			
			c=[0.2929 -0.3536 0.0884 -0.0442 0.0276 -0.0193];
			%coefficients in onsager-fuoss paper

		    factor=c*r';
			factor=factor';
		    mob_new=obj.F*omega-(obj.F*0.78420*z_list.*factor.*omega+31.410e-9).*sqrt(obj.I)./(1+1.5*sqrt(obj.I));
			mob_new=(mob_new.*z_list)';
			
			
		end
		
        
    end %End methods section
end %End class definition
