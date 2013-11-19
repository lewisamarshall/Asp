function [pass]=ion_tests(warning_flag)
	%ION_TESTS runs a set of ion test functions to test whether the ion class has bugs.
	if ~exist('warning_flag', 'var')
		warning off;
	end
	
	pass=[];
	pass=[pass, initialization_test()];
	pass=[pass, L_test()];
	pass=[pass, ionization_fraction_test()];
	
	if ~exist('warning_flag', 'var')
		warning on;
	end
end

function[pass]=initialization_test()
	try
		TRIS=ion('Tris', 1, 8.076, 29.5e-9);
		CL=ion('Chloride', -1, -2, -79.1e-9);
		ASCORBIC=ion('Ascorbic Acid', [-1, -2],[4.25, 10.5], [25.5e-9, 51e-9]);
		if ~all(ASCORBIC.absolute_mobility==(-[51e-9, 25.5e-9]))
			disp('Mobilities not corrected properly.')
			error('Wrong mobilities.')
		end
		
		pass=1;
		disp('Initialization test passed.')
	catch
		pass=0;
		disp('Initialization test failed.')
	end
end

function[pass]=L_test()
	try
		TEST=ion('TEST', [-2, -1, 1, 2], [-2, -1, 1, 2], [-2, -1, 1, 2]*1e-8);

		if ~all(TEST.z0==[-2, -1,0,  1, 2])
			disp('Problem inserting 0 in z0.')
			error('Problem inserting 0 in z0.')
		end

		if ~all(TEST.L==([1000, 10, 1, 10, 1000]))
			disp('L not calculated properly.')
			error('L not calculated properly.')
		end
		
		pass=1;
		disp('L test passed.')
	catch
		pass=0;
		disp('L test failed.')
	end
end

function [pass]=ionization_fraction_test()
	TEST1=ion('TEST ACID', -1, 7, -1);
	TEST2=ion('TEST BASE', 1, 7, 1);
	
	try
		if TEST1.ionization_fraction(7)~=0.5
			disp('Ionization fraction calculated incorrectly.')
			error('Ionization fraction calculated incorrectly.')
		end
		
		if TEST2.ionization_fraction(7)~=0.5
			disp('Ionization fraction calculated incorrectly.')
			error('Ionization fraction calculated incorrectly.')
		end
		
		TEST1.ionization_fraction(4, 0.1);
		TEST2.ionization_fraction(8, 0.1);
			
		pass=1;
		disp('Ionization test passed.')
	catch
		pass=0;
		disp('Ionization test failed.')
	end
end

function [pass]=Ka_test()
	pass=1;
end