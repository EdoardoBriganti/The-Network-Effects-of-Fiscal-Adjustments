function [data,T,n,tb,eb,dy_m,A,A_hat_tr,share_EB,share_TB] = data_upload(cutoff,num_ind)
%{
This function is meant to upload the data: 

INPUT      cutoff       scalar which must be equal to: a) 50; b) 70; c) 80
                        depending on whether you choose to construct EB
                        fiscal adjsutment shocks with at least 50%,70% or
                        80% of them coming from a spending cut. Default is
                        50%.
           industry     scalar which must be equal to: a) 15; b) 65; 
                        depending on the level of aggregation

OUTPUT     data         structure, which contains all the available data.
           T            scalar, time series length.
           n            scalar, cross section length.
           tb           vector (T x 1), is a dummy vector which identifies
                        in which year there is a TB fiscal adjustment plan
           eb           vector (T x 1), is a dummy vector which identifies
                        in which year there is an EB fiscal adjustment plan
           dy_m         matrix (n x T), contains the real growth rate of
                        value added of each sector. The first row
                        represents the time series for the first sector. 
           A            matrix (n x n), I-O matrix.
           A_hat_tr     matrix (n x n), I-O transformed matrix.
           wss          vector (n x 1), represent the vector of weights of
                        the spending shocks for each sector (since each
                        sector purchase from the government in a different
                        fashion.

%=========================================================================%

Favero, Karamysheva, Briganti 
"The network effects of fiscal adjustments" (2019)

%=========================================================================%
%}


%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Step 1) Upload Database:
        

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Database 1: Nominal Value Added from 1975 to 2014:

% Upload:
VA  = readmatrix('GDPbyInd_VA_1947-2015.xlsx','Sheet','VA','Range','C7:BS104');
VA_year = readmatrix('GDPbyInd_VA_1947-2015.xlsx','Sheet','VA','Range','C6:BS6');
VA_names = string(table2array(readtable('GDPbyInd_VA_1947-2015.xlsx','Sheet','VA','Range','B6:B104')));


switch num_ind 
    
    case 15
        
        % Update Industries and years:
        indices_ind = [find(VA_names=="Agriculture, forestry, fishing, and hunting");...
            find(VA_names=="Mining"); ...
            find(VA_names=="Utilities");...
            find(VA_names=="Construction");...
            find(VA_names=="Manufacturing");...
            find(VA_names=="Wholesale trade");...
            find(VA_names=="Retail trade");...
            find(VA_names=="Transportation and warehousing");...
            find(VA_names=="Information");...
            find(VA_names=="Finance, insurance, real estate, rental, and leasing");...
            find(VA_names=="Professional and business services");...
            find(VA_names=="Educational services, health care, and social assistance");...
            find(VA_names=="Arts, entertainment, recreation, accommodation, and food services");...
            find(VA_names=="Other services, except government");...
            find(VA_names=="Government enterprises")]; 
        
        
        
    case 19
        
        % Update Industries and years:
        indices_ind = [find(VA_names=="Agriculture, forestry, fishing, and hunting");...
            find(VA_names=="Mining"); ...
            find(VA_names=="Utilities");...
            find(VA_names=="Construction");...
            find(VA_names=="Manufacturing");...
            find(VA_names=="Wholesale trade");...
            find(VA_names=="Retail trade");...    
            find(VA_names=="Transportation and warehousing");...
            find(VA_names=="Warehousing and storage");...
            find(VA_names=="Information");...
            find(VA_names=="Finance and insurance");...
            find(VA_names=="Real estate and rental and leasing");...
            find(VA_names=="Professional and business services");...
            find(VA_names=="Educational services");...
            find(VA_names=="Health care and social assistance");...
            find(VA_names=="Arts, entertainment, and recreation");...
            find(VA_names=="Accommodation and food services");...
            find(VA_names=="Other services, except government");...
            find(VA_names=="Government enterprises")]; 
        
        % Subtract from "Transportation and Warehousing" the "Warehousing"
        % part.
        i = [find(VA_names=="Transportation and warehousing");...
            find(VA_names=="Warehousing and storage")];
        VA(i(1),:) = VA(i(1),:) - VA(i(2),:);
        % Rename Variable:
        VA_names(i(1)) = "Transportation";
        
        
    case 62
        
        % Update Industries and years:
        indices_ind = [find(VA_names=="Farms");...
            find(VA_names=="Forestry, fishing, and related activities");...
            find(VA_names=="Oil and gas extraction");...
            find(VA_names=="Mining, except oil and gas");...
            find(VA_names=="Support activities for mining");...
            find(VA_names=="Utilities");...
            find(VA_names=="Construction");...
            find(VA_names=="Wood products");...
            find(VA_names=="Nonmetallic mineral products");...
            find(VA_names=="Primary metals");...
            find(VA_names=="Fabricated metal products");...
            find(VA_names=="Machinery");...
            find(VA_names=="Computer and electronic products");...
            find(VA_names=="Electrical equipment, appliances, and components");...
            find(VA_names=="Motor vehicles, bodies and trailers, and parts");...
            find(VA_names=="Other transportation equipment");...
            find(VA_names=="Furniture and related products");...
            find(VA_names=="Miscellaneous manufacturing");...
            find(VA_names=="Food and beverage and tobacco products");...
            find(VA_names=="Textile mills and textile product mills");...
            find(VA_names=="Apparel and leather and allied products");...
            find(VA_names=="Paper products");...
            find(VA_names=="Printing and related support activities");...
            find(VA_names=="Petroleum and coal products");...
            find(VA_names=="Chemical products");...
            find(VA_names=="Plastics and rubber products");...
            find(VA_names=="Wholesale trade");...
            find(VA_names=="Retail trade");...
            find(VA_names=="Air transportation");...
            find(VA_names=="Rail transportation");...
            find(VA_names=="Water transportation");...
            find(VA_names=="Truck transportation");...
            find(VA_names=="Transit and ground passenger transportation");...
            find(VA_names=="Pipeline transportation");...
            find(VA_names=="Other transportation and support activities");...
            find(VA_names=="Warehousing and storage");...
            find(VA_names=="Publishing industries, except internet (includes software)");...
            find(VA_names=="Motion picture and sound recording industries");...
            find(VA_names=="Broadcasting and telecommunications");...
            find(VA_names=="Data processing, internet publishing, and other information services");...
            find(VA_names=="Federal Reserve banks, credit intermediation, and related activities");...
            find(VA_names=="Securities, commodity contracts, and investments");...
            find(VA_names=="Insurance carriers and related activities");...
            find(VA_names=="Funds, trusts, and other financial vehicles");...
            find(VA_names=="Real estate");...
            find(VA_names=="Rental and leasing services and lessors of intangible assets");...
            find(VA_names=="Legal services");...
            find(VA_names=="Computer systems design and related services");...
            find(VA_names=="Miscellaneous professional, scientific, and technical services");...
            find(VA_names=="Management of companies and enterprises");...
            find(VA_names=="Administrative and support services");...
            find(VA_names=="Waste management and remediation services");...
            find(VA_names=="Educational services");...
            find(VA_names=="Ambulatory health care services");...
            find(VA_names=="Hospitals and nursing and residential care facilities");...
            find(VA_names=="Social assistance");...
            find(VA_names=="Performing arts, spectator sports, museums, and related activities");...
            find(VA_names=="Amusements, gambling, and recreation industries");...
            find(VA_names=="Accommodation");...
            find(VA_names=="Food services and drinking places");...
            find(VA_names=="Other services, except government");...
            find(VA_names=="Government enterprises")];
        
end
% Update Industry Names: 
VA_names = VA_names(indices_ind);

% Indices of the years we are interested in: from 1975 to 2014:
indices_year = find(VA_year==1975):find(VA_year==2014);
% Update Years: 
VA_year = VA_year(:,indices_year);

% Update VA:
VA = VA(indices_ind,indices_year);


% Aggregate Government Enterprises:
% Identify rows to aggregate:
i  = find(VA_names=="Government enterprises");
% Aggregate by rows:
agg = sum(VA(i,:),1);
% Replace aggregated values in the original matrix in the same
% place of the first element to be aggregated:
VA(i(1),:) = agg;
% Eliminate the remaining rows:
i = i(2:end);
VA(i,:)     = [];
VA_names(i) = [];


% Make sure the desired number of sectors matches with the one constructed:
if size(VA,1)~=num_ind
    error('Aggregation Error')
end

% Data on Time Series:
data.years = VA_year;
T = length(data.years);    % Length of the time series




%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Database 2: 1997 I-O Matrices and Reliance to G of 1997

% Upload:
switch num_ind
    case 15
        load('IO_Database_15')
    case 19
        load('IO_Database_19')
    case 62
        load('IO_Database_62')
end
% Save names of industries:
data.industries = industries;

% Length of the cross section:
n = length(industries);  

% Make sure the number of sectors match:
if num_ind~=n
    error('Number of sectors mismatch')
end


% Construct Industry Shifter:
% 1) Reliance to Government as shifter for EB:
share_EB = Reliance2G;               % Dimension: (n x 1)
share_EB = share_EB ./ (ones(1,n)*share_EB);   % Normalization
% 2) Arithmetic Average for TB:
share_TB = 1/n .* ones(n,1);


        
        



%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Step 2) FISCAL PLANS
%  (Time series starts in year 1978 and ends in year 2014).

%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Database 3: Fiscal Consolidations.

% Upload:
Input = xlsread('CEM.xlsx','Input');

% Rename:
data.SUT = Input(1:end,2);
data.SAT = Input(1:end,3);
data.SAT1 = Input(1:end,4);
data.SAT2 = Input(1:end,5); 
data.SAT3 = Input(1:end,6);
data.SAT4 = Input(1:end, 7);
data.SAT5 = Input(1:end,8);
data.SUG = Input(1:end,9);
data.SAG = Input(1:end,10);
data.SAG1 = Input(1:end,11);
data.SAG2 = Input(1:end,12);
data.SAG3 = Input(1:end,13);
data.SAG4 = Input(1:end,14);
data.SAG5 = Input(1:end,15);

data.tax_shocks = [data.SUT data.SAT data.SAT1 data.SAT2 data.SAT3];
data.exp_shocks = [data.SUG data.SAG data.SAG1 data.SAG2 data.SAG3];


% Storing DUMMIES - First year is 1978 and last year is 2014
tb = Input(1:end,16);
eb_50 = Input(1:end,17);
eb_70 = Input(1:end,19);
eb_80 = Input(1:end,21);
switch cutoff
    case 50
        eb = eb_50;
    case 70
        eb = eb_70;
    case 80
        eb = eb_80;
end

% Computing FUTURE tax and spending shocks:
data.SAT_F = data.SAT1 + data.SAT2 + data.SAT3; 
data.SAG_F = data.SAG1 + data.SAG2 + data.SAG3;

% Computing total fiscal adjustments:
data.SU = data.SUT + data.SUG;       % Unexpected
data.SA = data.SAT + data.SAG;       % Announced
data.SA_F = data.SAT_F + data.SAG_F; % Future


 
%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Construction of Real Value Added:

% Upload Deflator:
data.deflator = readmatrix('CEM.xlsx','Sheet','VA','Range','B2:B41') ./ 100; % First year is 1975 - Last Year is 2014
% Rename Value Added:
data.y_nominal = VA;      % First year is 1975
% Construct Real Value Added:
data.y_real  = data.y_nominal./transpose(data.deflator);    % First year is 1975


% Constructing percentage change of sectors real value added
data.dy_real = ( diff(data.y_real,1,2) ./ data.y_real(:,1:end-1) ) .* 100; % First year/column is 1976 
dy_m     = data.dy_real(:,3:end);   % matrix (n x T) - first year 1978, last year 2014

% Aggregate data:
data.y_real_USA  = sum(data.y_real,2); % this is a (40 x 1) vector (first year 1975)
data.dy_real_USA = diff(data.y_real_USA) ./ data.y_real_USA(1:end-1) * 100; % first year is 1976
% Since we are interested only in years from 1977 to 2014 (because the
% first shock is in year 1978 and the model is lagged - therefore we need
% the value for 1977)
% We can drop the first observation:
data.dy_real_USA = data.dy_real_USA(2:end); % first year 1977


% Data on Time Series:
data.years = data.years(find(data.years==1978):find(data.years==2014));
T = length(data.years);    % Length of the time series


%-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Construct Industry Shares: from 1975 to 2014
y_share = data.y_nominal ./ sum(data.y_nominal,1); % matrix (n x 40)
% Select shares only in the years from 1978 to 2014)
y_share = y_share(:,4:end) ;
% Take average of industry share conditional on the occurrence of a shock:
% TB:
ind_weight_TB = y_share .* transpose(tb); 
ind_weight_TB(:,~any(ind_weight_TB)) = [];
data.ind_weight_TB = mean(ind_weight_TB,2); % (n x 1)
% EB:
ind_weight_EB = y_share .* transpose(eb); 
ind_weight_EB(:,~any(ind_weight_EB)) = [];
data.ind_weight_EB = mean(ind_weight_EB,2); % (n x 1)

end 