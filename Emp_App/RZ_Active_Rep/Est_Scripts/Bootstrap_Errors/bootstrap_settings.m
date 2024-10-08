% settings for bootstrapped standard errors

% DWB settings
e1 = 1; 
e2 = -1; 
p1 = 0.5; 
p2 = 1-p1;
values = [e1, e2]; 
probabilities = [p1, p2]; 
dwb_settings = [values; probabilities];

bw = floor(505^(1/4));