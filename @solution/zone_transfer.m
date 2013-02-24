function Qi=zone_transfer(obj, vol)
    %Calculates the zone transfer charge of the solution at a given
    %volume. The volume of the solution is specified in liters. 
    Qi=obj.concentrations./obj.get_transference.*vol.*obj.F;
end