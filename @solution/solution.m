classdef solution < handle
    % A for a solution containing one or more ions.
	% Upon initialization, immediately determines the pH and ionic strength of the solution. 
	% Contains methods for simulation important solution properties. 
    
    properties(Constant = true, GetAccess = 'private')
        F=9.65E4;           % Faraday's const.[C/mol]
        Rmu=8.31;           % Universal gas const. [J/mol*K]
        Temp=298;           % Temperature [K]
        Kw=1E-14;           % Water equilibrium constant
        muH=362E-9%/obj.F;   % Mobility of Hydronium   % [m^2/s*V]/F --> [mol*s/Kg]
        muOH=205E-9%/obj.F;  % Mobility of Hydroxide   % [m^2/s*V]/F --> [mol*s/Kg]
        Lpm3=1000;          % Liters per meter^3
        visc=1E-3;          % Dynamic viscosity (water) [Pa s]
		Adh=0.512; 			% L^1/2 / mol^1/2, approximate for room temperature
		aD=1.5*sqrt(2); 	% mol^-1/2 mol^-3/2, approximation
    end
    
    properties
        ions;
        concentrations;
        pH;
        % cond; % conductivity is no longer a static object property.
		I=0;
    end
    
    methods
        
        function obj=solution(ions, concentrations)
			%Class Constructor

            if(nargin == 2)
				if ~iscell(ions)
					ions=num2cell(ions);
				end
				
                if isvector(ions) && all(strcmp(cellfun(@class, ions, 'UniformOutput', false), 'ion'))
                    obj.ions=ions;
                else
                    error('You must input a cell vector of ions.')
                end
                if isvector(concentrations) && length(ions)==length(concentrations)
					if iscell(concentrations)
						concentrations=cell2mat(concentrations)
					end
					obj.concentrations=concentrations;
                else
                    error('The concentrations vector must be the same size as the ions vector.')
                end
                
            else
                error('Solutions must have a cell of ions and a cell or vector of concentrations.')
            end
            
            obj.pH=get_pH(obj);
            % obj.cond=get_conductivity(obj);
			obj.I=ionic_strength(obj, obj.pH);
        end

		function new_sol=add_ion(obj, new_ion, new_conc)
			new_sol=solution(cat(2, obj.ions, {new_ion}), cat(2, obj.concentrations, new_conc));
		end
        
        function pH=get_pH(obj, I)
			% Finds the pH of the object. If an ionic strength is specified, 
			% uses the corrected acidity constants. 
			
			% If ionic strength does not exist, assume it is zero.
			% This function is used to find the equilibrium state, 
			% so it cannot pull the ionic strength from the object.
			if ~exist('I', 'var')
				I=0; 
			end
			
			%Set up concentration vector
            cMat=obj.concentrations;
			
			% Find the order of the polynomial. This is the maximum
			% size of the list of charge states in an ion. 
            MaxCol=-inf;
            for i=1:length(obj.concentrations)
                MaxCol=max(MaxCol, max(obj.ions{i}.z)-min(obj.ions{i}.z)+1);
            end
            
			% Set up the matrix of Ls, the multiplication
			% of acidity coefficients for each ion.
            LMat=zeros(length(cMat),MaxCol);
            
            for i=1:length(cMat)
                L_list=obj.ions{i}.get_L(I);
                Ip1=find(obj.ions{i}.z==1);     Im1=find(obj.ions{i}.z==-1);
                L_list=[L_list(1:Im1),1,L_list(Ip1:end)];
                LMat(i,1:length(obj.ions{i}.z)+1)=L_list;
            end
            
        
            % Construct Q 
			% Made more consice by removing Q1 and Q2.  
			Q=1;
			% Convolve each line of the L matrix. 
            for j=1:size(LMat,1)
                Q=conv(Q,LMat(j,:));
            end %for j
			%Convolve with water dissociation.
			Q=conv(Q, [-obj.Kw_eff(I) 0 1]);
			
			
            
            
            %%%
            %Construct P (done)
            %%%
%             PMat=zeros(obj.concentrations, 1)
            for i=1:length(obj.concentrations)
                z_list=obj.ions{i}.get_z0;
                
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
            
            
            %%%
            %Pad Smaller Matrix
            %%%
            
            SizeDiff=size(Q,2)-size(PMat,2);
            if SizeDiff>0
                PMat=[PMat,repmat(PMat(:,1)*0,1,SizeDiff)];
            elseif SizeDiff<0
                Q=[Q,repmat(0,1,SizeDiff)];
            end
                 
         
            %%%
            %Construct Polynomial
            %%%
            cTotRep=cMat'*ones(1,size(PMat,2));

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
		
		function [pH, I]=find_equilibrium(obj)
			% Uses the fzero root finder to find the equilibrium
			% pH and ionic strength of a solution, using the ionic-
			% strength-adjusted activity coefficients.
			
			I=obj.ionic_strength(obj.get_pH);
			I=fzero(@(x)equil_offset(obj, x), I);
			pH=obj.get_pH(I);
		end
        
		function Residual=equil_offset(obj,I_i)
			% This function finds the offset between proposed
			% pH and I and the pH and I calculated from the initial guess.
			pH=obj.get_pH(I_i);
			I_f=obj.ionic_strength(pH, I_i);
			
			Residual=(I_f-I_i);
		end
		
        function cond=cond(obj)
            %Calculates the conductivity of the solution. Does not
            %currently include the conductivity contribution of protons or
            %hydroxyls. 
            cond=0;
            for i=1:length(obj.concentrations)
                cond=cond+obj.ions{i}.molar_conductivity(obj.pH)*obj.concentrations(i);
            end
            cond=cond+obj.F*obj.muH*obj.cH*obj.Lpm3; % H+ conductivity
            cond=cond+obj.F*obj.muOH*obj.cOH*obj.Lpm3; % H+ conductivity
        end
        
        function T=get_transference(obj)
            %Gets the fraction of charge carried by each of the ions.
            %Should not precisely add to 1, because some charge is carried
            %by protons and hydroxyls. 
            T=zeros(1,length(obj.ions));
            for i=1:length(T)
                T(i)=obj.ions{i}.molar_conductivity(obj.pH).*obj.concentrations(i);
            end
            T=T/obj.cond;
        end
        
        function Qi=passage_charge(obj, LE)
            Qi=0;
        end
        
        function I=ionic_strength(obj, pH, I_guess)
            % Calculates the ionic strength of the solution. This function
            % does not account for changes in activity due to ionic strength
            % effects. 
			
			% If a pH is not specified in the function call, use the pH
			% of the object.
			if ~exist('pH', 'var')
				pH=obj.pH;
			end
			
			if ~exist('I_guess', 'var') %calculate without an initial guess for ionic strength
            	I=0;
            	for i=1:length(obj.ions);
            		I=I+obj.concentrations(i)*sum((double(obj.ions{i}.z).^2).*obj.ions{i}.ionization_fraction(pH));
            	end
				% Add the ionic strength due to water dissociation. 
				I=I+obj.cH + obj.cOH;
            	I=I/2;
			else
            	I=0;
            	for i=1:length(obj.ions);
            		I=I+obj.concentrations(i)*sum((double(obj.ions{i}.z).^2).*obj.ions{i}.ionization_fraction(obj.pH, I_guess));
            	end
				I=I+obj.cH + obj.cOH;
            	I=I/2;
			end
        end
		        
        function Qi=zone_transfer(obj, vol)
            %Calculates the zone transfer charge of the solution at a given
            %volume. The volume of the solution is specified in liters. 
            Qi=obj.concentrations./obj.get_transference.*vol.*obj.F;
        end

		function Cb=buffering_capacity(obj)
			c=obj.concentrations(obj.concentrations>0);
			c=0.01*min(c);
			new_sol=obj.add_ion(ion('A-', -1, -2, -1), c);
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
			mu=conc_list.*z_list.^2./obj.ionic_strength/2; %total potential
			
			
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
		    mob_new=obj.F*omega-(obj.F*0.78420*z_list.*factor.*omega+31.410e-9).*sqrt(obj.ionic_strength)./(1+1.5*sqrt(obj.ionic_strength));
			mob_new=(mob_new.*z_list)';
			
			
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
		
		function cH=cH(obj)
			% Supplies the concentration of H+ in the solution. 
			cH=10^(-obj.pH);
		end
		
		function cOH=cOH(obj)
			% Supplies the concentration of OH- in the solution. 
			cOH=obj.Kw_eff/obj.cH;
		end
			
        
    end %End methods section
end %End class definition

% Outdated code.

% Q1=1;
% for j=1:size(LMat,1)
%     Q1=conv(Q1,LMat(j,:));
% end %for j
% Q2=[-obj.Kw_eff(I) 0 1];
% Q=conv(Q1,Q2);