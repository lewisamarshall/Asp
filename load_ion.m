function [loaded_ion]=load_ion(ion_name)
	% LOAD_ION is a function to pull an ion by name from a database.
	% ion=LOAD_ION('ion name') pulls the named ion from the database.
	% Database derived from Peakmaster.
	
	load database.mat
	N=find(strcmpi(ion_name, NAME));
	
	if isempty(N)
		error('Ion not found in database.')
	end
	
	indices=~isnan(PKA(N, :));
	loaded_ion=ion( NAME{N},...
		Z(indices),...
		PKA(N, indices),...
		MOBILITY(N, indices).*1e-9.*sign(Z(indices)));
end