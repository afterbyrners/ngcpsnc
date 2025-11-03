function CL_a = polhamus(AR,lambda,sweep,M) % CL_alpha function 
   if AR < 4
       k = 1 + (AR*(1.87-0.000233*sweep))/100;
   elseif AR >= 4
       k = 1+((8.2 -2.3*sweep)-AR*(0.22-0.153*sweep))/100;
   end
   tangent_half = tan(sweep)-(4*0.5*(1-lambda))/(AR*(1+lambda));
   CL_a = 2*pi*AR/(2+sqrt(((AR^2)*(1-M^2))/(k^2)*(1 + ((tangent_half)^2)/(1-M^2))+4));
end