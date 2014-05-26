function [p]=separability(ion1, ion2, pH, I)
	% Separability is a function to calculate the separability between two
	% ions. The second ion is used in the denominator. The function returns
	% a signed quantity, but the unsigned separability can be found by using abs(separability(...)).
	
	if strcmp(class(pH), 'solution')
		SOL=pH;
		pH=SOL.pH;
		I=SOL.I;
	end
	if ~exist('I', 'var')
		I=0;
	end
	
	p=(ion1.effective_mobility(pH, I)-ion2.effective_mobility(pH, I))/ion2.effective_mobility(pH, I);
	
end