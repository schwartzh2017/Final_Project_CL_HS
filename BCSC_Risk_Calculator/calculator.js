// Copyright 2012-2015, Breast Cancer Surveillance Consortium
// Developed by Jambeyang Research, LLC., http://www.jambeyang.com, email: wlin@jambeyang.com
// Last modified on 5/24/2015

function include(file)
{
  var script  = document.createElement('script');
  script.src  = file;
  script.type = 'text/javascript';
  script.defer = true;
  document.getElementsByTagName('head').item(0).appendChild(script);
}
include('mapBiopsy.js');
include('mapDensity.js');
include('mapFamhx.js');

var risk_5Y = 0.0;
var risk_10Y = 0.0;
var risk_5Y_average = 0.0;
var risk_10Y_average = 0.0;
var risk_avg_low = 0.0;
var risk_avg_high = 0.0;
var risk_top5 = 0.0;

var age_var;
var e = 2.718281828;

var low_cutoff = 0.75;
var high_cutoff = 1.25;

// race = 1 for white, 2 for black, 3 for asian, 4 for native_american, 5 for hispanic, 6 for unknown_and_other
// relativeBC = 0 for No, 1 for Yes
// biopsy = 1 for unknown, 2 for No, 3 for prior with unknown diagnosis, 4 for unproliferative,
//          5 for proliferative w/o atypia, 6 for atypia, 7 for LCIS
function calculate_risk (age, race, relativeBC, biopsy, density) {
	//alert("age=" + age + " race=" + race + " relativeBC=" + relativeBC + " biopsy=" + biopsy + " density=" + density);
	var temprace = race;	// OK for w-white, 2-black, 3-asian
	if (race == 1)
	{
		temprace = 1; 	// white
	}
	else if (race == 2)
	{
		temprace = 2; 	// black
	}
	else if (race == 3)
	{
		temprace = 3; 	// asian
	} 
	else if (race == 4)
	{
		temprace = 5; 	// native american *** change ***
	} 
	else if (race == 5)
	{
		temprace = 8;	// hispanic *** change ***
	}
	else
	{
		temprace = 1;	// for rest, treat as white
	}
	computerRiskV2(age, temprace, relativeBC, biopsy, density);	// version 2
	risk_5Y = risk_5Y.toFixed(2);
	compute_average_risk_by_age_race (age, race);	// version 2
}

// Compute risk model 2.0
// preliminary version 5/14/2015
function computerRiskV2(age, temprace, fam_history, biopsy_result, birads_density)
{
	//alert("age=" + age + " temprace=" + temprace + " fam_history=" + fam_history + " biopsy_result=" + biopsy_result + " birads_density=" + birads_density);	// for debug
	var adjustment_density_arg = (1000 * temprace) + (10 * age) + birads_density;
	//alert (adjustment_density_arg);
	var adjustment_density = map_adjustment_density(adjustment_density_arg);	// defined in mapDensity.js
	//alert (adjustment_density);	// for debug
	var adjustment_famhx_arg = (1000 * temprace) + (10 * age) + fam_history;
	var adjustment_famhx = map_adjustment_famhx(adjustment_famhx_arg);			// defined in mapFamhx.js
	//alert (adjustment_famhx);		// for debug
	var adjustment_biopsy_arg = (1000 * temprace) + (10 * age) + biopsy_result;
	var adjustment_biopsy = map_adjustment_biopsy(adjustment_biopsy_arg);		// defined in mapBiopsy.js
	//alert (adjustment_biopsy);		// for debug
	
	// initialize time fields
	var inc_times = new Array(12); 
    for (var iloop = 0; iloop < 12; iloop++)
    {
        inc_times[iloop] = 0;
    }
	
	// initialize variables;
	var prop_notatrisk = 0;
	var prop_atrisk = 1;
	
	// for 5-year and 10-year risk calculation
	var loopend = 10;
	
	// iteration 1
	for (var i = 1; i <= loopend; i++)
    {
        var agemod = age + i - 1;
      
        if ((temprace == 1) || (temprace == 6))
        {
            // WHITE, NON HISPANIC OR OTHER, MIXED (2+ RACES), UNKNOWN
            inc_times[i] = adjustment_density * adjustment_famhx * adjustment_biopsy *
                            ((-0.000007162 * (agemod * agemod * agemod)) 
                            + (0.001111136 * (agemod * agemod)) 
                            - (0.043528999 * agemod) 
                            + 0.514780021);
        }
        else if (temprace == 2)
        {
            // BLACK, NON HISPANIC
            inc_times[i] = adjustment_density * adjustment_famhx * adjustment_biopsy * 
                            ((-0.000004332 * (agemod * agemod * agemod)) 
                            + (0.000644475 * (agemod * agemod)) 
                            - (0.021243209 * agemod) 
                            + 0.195808387);
        }
        else if (temprace == 3)
        {
            // ASIAN, NATIVE HAWAIIAN, OR PACIFIC ISLANDER
            inc_times[i] = adjustment_density * adjustment_famhx * adjustment_biopsy * 
                            ((-0.000003002 * (agemod * agemod * agemod)) 
                            + (0.000372751 * (agemod * agemod)) 
                            - (0.007182055 * agemod) 
                            - 0.025628846);
        }
        else if (temprace == 5)
        {
            // AMERICAN INDIAN OR ALASKA NATIVE
            inc_times[i] = adjustment_density * adjustment_famhx * adjustment_biopsy * 
                            ((-0.000005880 * (agemod * agemod * agemod)) 
                            + (0.000861287 * (agemod * agemod)) 
                            - (0.032530311 * agemod) 
                            + 0.373493536);
        }
        else if (temprace == 8)
        {
            // HISPANIC
            inc_times[i] = adjustment_density * adjustment_famhx * adjustment_biopsy * 
                            ((-0.000004211 * (agemod * agemod * agemod)) 
                            + (0.000628530 * (agemod * agemod)) 
                            - (0.022484578 * agemod) 
                            + 0.234323344);
        }
    }

	// iteration 2
	for (var j = 2; j <= loopend; j++)
    {
        var agemod = age + j - 2;
        if ((temprace == 1) || (temprace == 6))
        {
            // WHITE, NON HISPANIC OR OTHER, MIXED (2+ RACES), UNKNOWN
            prop_notatrisk = prop_notatrisk 
                            + prop_atrisk * (0.004741239 * Math.pow(e, 0.082473361 * agemod)) / 100
                            + prop_atrisk * adjustment_density * adjustment_famhx * adjustment_biopsy 
                                * ((-0.000002094 * (agemod * agemod * agemod)) 
                                + (0.000284086 * (agemod * agemod)) 
                                - (0.009292718 * agemod) + 0.080659188) 
                                / 100
                            + prop_atrisk * adjustment_density * adjustment_famhx * adjustment_biopsy
                                * ((-0.000007162 * (agemod * agemod * agemod)) 
                                + (0.001111136 * (agemod * agemod)) 
                                - (0.043528999 * agemod) + 0.514780021) 
                                / 100;
        }
        else if (temprace == 2)
        {
            // BLACK, NON HISPANIC
            prop_notatrisk = prop_notatrisk
                            + prop_atrisk * (0.008976903 * Math.pow(e, 0.077540550 * agemod)) / 100
                            + prop_atrisk * adjustment_density * adjustment_famhx * adjustment_biopsy
                                * ((-0.000002706 * (agemod * agemod * agemod))
                                + (0.000403385 * (agemod * agemod)) 
                                - (0.016548963 * agemod) + 0.207686538) 
                                / 100
                            + prop_atrisk * adjustment_density * adjustment_famhx * adjustment_biopsy
                                * ((-0.000004332 * (agemod * agemod * agemod))
                                + (0.000644475 * (agemod * agemod)) 
                                - (0.021243209 * agemod) + 0.195808387) 
                                / 100;                        
        }
        else if (temprace == 3)
        {
            // ASIAN, NATIVE HAWAIIAN, OR PACIFIC ISLANDER
            prop_notatrisk = prop_notatrisk + prop_atrisk * (0.001324086 * Math.pow(e, 0.091167044 * agemod)) / 100
                            + prop_atrisk * adjustment_density * adjustment_famhx * adjustment_biopsy
                                * ((-0.000001496 * (agemod * agemod * agemod))
                                + (0.000185804 * (agemod * agemod)) 
                                - (0.004747568 * agemod) + 0.016457265) 
                                / 100
                            + prop_atrisk * adjustment_density * adjustment_famhx * adjustment_biopsy
                                * ((-0.000003002 * (agemod * agemod * agemod))
                                + (0.000372751 * (agemod * agemod)) 
                                - (0.007182055 * agemod) - 0.025628846) 
                                / 100;
        }
        else if (temprace == 5)
        {
            // AMERICAN INDIAN OR ALASKA NATIVE
            prop_notatrisk = prop_notatrisk + prop_atrisk * (0.016491843 * Math.pow(e, 0.067539660 * (agemod))) / 100
                            + prop_atrisk * adjustment_density * adjustment_famhx * adjustment_biopsy
                                * ((0.000000191 * (agemod * agemod * agemod))
                                - (0.000049665 * (agemod * agemod)) 
                                + (0.004730528 * agemod) - 0.102520406) 
                                / 100
                            + prop_atrisk * adjustment_density * adjustment_famhx * adjustment_biopsy
                                * ((-0.000005880 * (agemod * agemod * agemod))
                                + (0.000861287 * (agemod * agemod)) 
                                - (0.032530311 * agemod) + 0.373493536) 
                                / 100;
        }
        else if (temprace == 8)
        {
            // HISPANIC
            prop_notatrisk = prop_notatrisk + prop_atrisk * (0.002253560 * Math.pow(e, 0.088334802 * (agemod))) / 100
                            + prop_atrisk * adjustment_density * adjustment_famhx * adjustment_biopsy
                                * ((-0.000001403 * (agemod * agemod * agemod))
                                + (0.000195157 * (agemod * agemod)) 
                                - (0.006832187 * agemod) + 0.067237073) 
                                / 100
                            + prop_atrisk * adjustment_density * adjustment_famhx * adjustment_biopsy
                                * ((-0.000004211 * (agemod * agemod * agemod))
                                + (0.000628530 * (agemod * agemod)) 
                                - (0.022484578 * agemod) + 0.234323344) 
                                / 100;
        }

        prop_atrisk = 1 - prop_notatrisk;
        inc_times[j] = inc_times[j] * prop_atrisk;
    }
	
	// compute 5 year risk
    var agerace_inc_time_5 = 0;
    for (var i = 1; i <= 5; i++)
    {
        agerace_inc_time_5 += inc_times[i];
    }

    // compute 10 year risk
    var agerace_inc_time_10 = 0;
    for (var i = 1; i <= 10; i++)
    {
        agerace_inc_time_10 += inc_times[i];
    }
	
	// save result
	risk_5Y = agerace_inc_time_5;
	risk_10Y = agerace_inc_time_10;
}


function compute_average_risk_by_age_race (age, race)
{
	risk_5Y_average = 0.0;
	risk_10Y_average = 0.0;
	if (age < 35) 
	{
		risk_5Y_average = 0.0;
	}
	else if (age > 74)
	{
		risk_5Y_average = 0.0;
	}
	else
	{
		if (race == 2)	// black
		{
			if (age == 35)
			{
				risk_5Y_average = 0.3631726297;
				risk_10Y_average = 0.9499765133;
			}
			else if (age == 36)
			{
				risk_5Y_average = 0.4067185163;
				risk_10Y_average = 1.0415708561;
			}
			else if (age == 37)
			{
				risk_5Y_average = 0.45170812;
				risk_10Y_average = 1.1351932969;
			}
			else if (age == 38)
			{
				risk_5Y_average = 0.4980057896;
				risk_10Y_average = 1.2305584057;
			}
			else if (age == 39)
			{
				risk_5Y_average = 0.5454757079;
				risk_10Y_average = 1.3273801884;
			}
			else if (age == 40)
			{
				risk_5Y_average = 0.5939818912;
				risk_10Y_average = 1.4253720617;
			}
			else if (age == 41)
			{
				risk_5Y_average = 0.6433881796;
				risk_10Y_average = 1.5242468174;
			}
			else if (age == 42)
			{
				risk_5Y_average = 0.6935582275;
				risk_10Y_average = 1.6237165791;
			}
			else if (age == 43)
			{
				risk_5Y_average = 0.7443554904;
				risk_10Y_average = 1.723492751;
			}
			else if (age == 44)
			{
				risk_5Y_average = 0.7956432079;
				risk_10Y_average = 1.823285963;
			}
			else if (age == 45)
			{
				risk_5Y_average = 0.847284386;
				risk_10Y_average = 1.9228060112;
			}
			else if (age == 46)
			{
				risk_5Y_average = 0.8991417756;
				risk_10Y_average = 2.0217617995;
			}
			else if (age == 47)
			{
				risk_5Y_average = 0.95107785;
				risk_10Y_average = 2.1198612831;
			}
			else if (age == 48)
			{
				risk_5Y_average = 1.0029547804;
				risk_10Y_average = 2.2168114182;
			}
			else if (age == 49)
			{
				risk_5Y_average = 1.0546344104;
				risk_10Y_average = 2.3123181215;
			}
			else if (age == 50)
			{
				risk_5Y_average = 1.1059782307;
				risk_10Y_average = 2.4060862451;
			}
			else if (age == 51)
			{
				risk_5Y_average = 1.156847353;
				risk_10Y_average = 2.4978195702;
			}
			else if (age == 52)
			{
				risk_5Y_average = 1.2071024857;
				risk_10Y_average = 2.5872208276;
			}
			else if (age == 53)
			{
				risk_5Y_average = 1.2566039115;
				risk_10Y_average = 2.6739917501;
			}
			else if (age == 54)
			{
				risk_5Y_average = 1.3052114678;
				risk_10Y_average = 2.7578331657;
			}
			else if (age == 55)
			{
				risk_5Y_average = 1.3527845314;
				risk_10Y_average = 2.8384451398;
			}
			else if (age == 56)
			{
				risk_5Y_average = 1.3991820094;
				risk_10Y_average = 2.9155271749;
			}
			else if (age == 57)
			{
				risk_5Y_average = 1.4442623372;
				risk_10Y_average = 2.9887784818;
			}
			else if (age == 58)
			{
				risk_5Y_average = 1.4878834855;
				risk_10Y_average = 3.0578983296;
			}
			else if (age == 59)
			{
				risk_5Y_average = 1.5299029798;
				risk_10Y_average = 3.1225864926;
			}
			else if (age == 60)
			{
				risk_5Y_average = 1.5701779328;
				risk_10Y_average = 3.1825438042;
			}
			else if (age == 61)
			{
				risk_5Y_average = 1.6085650945;
				risk_10Y_average = 3.2374728372;
			}
			else if (age == 62)
			{
				risk_5Y_average = 1.6449209216;
				risk_10Y_average = 3.2870787241;
			}
			else if (age == 63)
			{
				risk_5Y_average = 1.6791016709;
				risk_10Y_average = 3.3310701385;
			}
			else if (age == 64)
			{
				risk_5Y_average = 1.7109635199;
				risk_10Y_average = 3.3691604541;
			}
			else if (age == 65)
			{
				risk_5Y_average = 1.7403627202;
				risk_10Y_average = 3.4010691035;
			}
			else if (age == 66)
			{
				risk_5Y_average = 1.7671557873;
				risk_10Y_average = 3.4265231559;
			}
			else if (age == 67)
			{
				risk_5Y_average = 1.791199733;
				risk_10Y_average = 3.4452591338;
			}
			else if (age == 68)
			{
				risk_5Y_average = 1.8123523474;
				risk_10Y_average = 3.4570250909;
			}
			else if (age == 69)
			{
				risk_5Y_average = 1.8304725361;
				risk_10Y_average = 3.4615829675;
			}
			else if (age == 70)
			{
				risk_5Y_average = 1.8454207199;
				risk_10Y_average = 3.4587112423;
			}
			else if (age == 71)
			{
				risk_5Y_average = 1.8570593065;
				risk_10Y_average = 3.4482078931;
			}
			else if (age == 72)
			{
				risk_5Y_average = 1.8652532416;
				risk_10Y_average = 3.4298936751;
			}
			else if (age == 73)
			{
				risk_5Y_average = 1.8698706493;
				risk_10Y_average = 3.4036157209;
			}
			else if (age == 74)
			{
				risk_5Y_average = 1.8707835715;
				risk_10Y_average = 3.3692514552;
			}
		}
		else if (race ==3) // Asian
		{
			if (age == 35)
			{
				risk_5Y_average = 0.3339135888;
				risk_10Y_average = 0.8675621631;
			}
			else if (age == 36)
			{
				risk_5Y_average = 0.3741590225;
				risk_10Y_average = 0.9478540205;
			}
			else if (age == 37)
			{
				risk_5Y_average = 0.4146819437;
				risk_10Y_average = 1.0281744382;
			}
			else if (age == 38)
			{
				risk_5Y_average = 0.4553912757;
				risk_10Y_average = 1.108339793;
			}
			else if (age == 39)
			{
				risk_5Y_average = 0.4961959995;
				risk_10Y_average = 1.1881666804;
			}
			else if (age == 40)
			{
				risk_5Y_average = 0.5370051464;
				risk_10Y_average = 1.2674718727;
			}
			else if (age == 41)
			{
				risk_5Y_average = 0.5777277899;
				risk_10Y_average = 1.3460722718;
			}
			else if (age == 42)
			{
				risk_5Y_average = 0.6182730348;
				risk_10Y_average = 1.4237848576;
			}
			else if (age == 43)
			{
				risk_5Y_average = 0.6585500067;
				risk_10Y_average = 1.500426632;
			}
			else if (age == 44)
			{
				risk_5Y_average = 0.6984678385;
				risk_10Y_average = 1.5758145597;
			}
			else if (age == 45)
			{
				risk_5Y_average = 0.7379356574;
				risk_10Y_average = 1.6497655061;
			}
			else if (age == 46)
			{
				risk_5Y_average = 0.776862569;
				risk_10Y_average = 1.7220961732;
			}
			else if (age == 47)
			{
				risk_5Y_average = 0.8151576424;
				risk_10Y_average = 1.792623036;
			}
			else if (age == 48)
			{
				risk_5Y_average = 0.8527298928;
				risk_10Y_average = 1.8611622778;
			}
			else if (age == 49)
			{
				risk_5Y_average = 0.8894882646;
				risk_10Y_average = 1.9275297299;
			}
			else if (age == 50)
			{
				risk_5Y_average = 0.9253416139;
				risk_10Y_average = 1.9915408146;
			}
			else if (age == 51)
			{
				risk_5Y_average = 0.9601986906;
				risk_10Y_average = 2.0530104951;
			}
			else if (age == 52)
			{
				risk_5Y_average = 0.9939681216;
				risk_10Y_average = 2.1117532349;
			}
			else if (age == 53)
			{
				risk_5Y_average = 1.026558394;
				risk_10Y_average = 2.1675829701;
			}
			else if (age == 54)
			{
				risk_5Y_average = 1.0578778405;
				risk_10Y_average = 2.2203130977;
			}
			else if (age == 55)
			{
				risk_5Y_average = 1.0878346256;
				risk_10Y_average = 2.2697564849;
			}
			else if (age == 56)
			{
				risk_5Y_average = 1.1163367358;
				risk_10Y_average = 2.3157255044;
			}
			else if (age == 57)
			{
				risk_5Y_average = 1.1432919722;
				risk_10Y_average = 2.3580321009;
			}
			else if (age == 58)
			{
				risk_5Y_average = 1.1686079485;
				risk_10Y_average = 2.3964878962;
			}
			else if (age == 59)
			{
				risk_5Y_average = 1.1921920939;
				risk_10Y_average = 2.4309043407;
			}
			else if (age == 60)
			{
				risk_5Y_average = 1.2139516645;
				risk_10Y_average = 2.4610929191;
			}
			else if (age == 61)
			{
				risk_5Y_average = 1.2337937616;
				risk_10Y_average = 2.4868654218;
			}
			else if (age == 62)
			{
				risk_5Y_average = 1.2516253622;
				risk_10Y_average = 2.5080342918;
			}
			else if (age == 63)
			{
				risk_5Y_average = 1.2673533614;
				risk_10Y_average = 2.5244130618;
			}
			else if (age == 64)
			{
				risk_5Y_average = 1.2808846303;
				risk_10Y_average = 2.5358168935;
			}
			else if (age == 65)
			{
				risk_5Y_average = 1.2921260916;
				risk_10Y_average = 2.5420632373;
			}
			else if (age == 66)
			{
				risk_5Y_average = 1.3009848173;
				risk_10Y_average = 2.5429726287;
			}
			else if (age == 67)
			{
				risk_5Y_average = 1.3073681503;
				risk_10Y_average = 2.5383696419;
			}
			else if (age == 68)
			{
				risk_5Y_average = 1.3111838563;
				risk_10Y_average = 2.5280840213;
			}
			else if (age == 69)
			{
				risk_5Y_average = 1.312340309;
				risk_10Y_average = 2.5119520135;
			}
			else if (age == 70)
			{
				risk_5Y_average = 1.3107467156;
				risk_10Y_average = 2.4898179257;
			}
			else if (age == 71)
			{
				risk_5Y_average = 1.306313387;
				risk_10Y_average = 2.4615359345;
			}
			else if (age == 72)
			{
				risk_5Y_average = 1.2989520612;
				risk_10Y_average = 2.4269721729;
			}
			else if (age == 73)
			{
				risk_5Y_average = 1.2885762867;
				risk_10Y_average = 2.3860071219;
			}
			else if (age == 74)
			{
				risk_5Y_average = 1.2751018753;
				risk_10Y_average = 2.3385383321;
			}
		}
		else if (race == 4) // Native american
		{
			if (age == 35)
			{
				risk_5Y_average = 0.2563190665;
				risk_10Y_average = 0.7023183227;
			}
			else if (age == 36)
			{
				risk_5Y_average = 0.2921043626;
				risk_10Y_average = 0.7802241028;
			}
			else if (age == 37)
			{
				risk_5Y_average = 0.3297468905;
				risk_10Y_average = 0.8607913932;
			}
			else if (age == 38)
			{
				risk_5Y_average = 0.3690657245;
				risk_10Y_average = 0.9436485876;
			}
			else if (age == 39)
			{
				risk_5Y_average = 0.4098799011;
				risk_10Y_average = 1.028424263;
			}
			else if (age == 40)
			{
				risk_5Y_average = 0.4520084359;
				risk_10Y_average = 1.1147472583;
			}
			else if (age == 41)
			{
				risk_5Y_average = 0.4952703376;
				risk_10Y_average = 1.2022467517;
			}
			else if (age == 42)
			{
				risk_5Y_average = 0.5394846227;
				risk_10Y_average = 1.2905523407;
			}
			else if (age == 43)
			{
				risk_5Y_average = 0.5844703287;
				risk_10Y_average = 1.3792941248;
			}
			else if (age == 44)
			{
				risk_5Y_average = 0.630046528;
				risk_10Y_average = 1.4681027954;
			}
			else if (age == 45)
			{
				risk_5Y_average = 0.6760323413;
				risk_10Y_average = 1.5566097334;
			}
			else if (age == 46)
			{
				risk_5Y_average = 0.7222469526;
				risk_10Y_average = 1.6444471196;
			}
			else if (age == 47)
			{
				risk_5Y_average = 0.768509625;
				risk_10Y_average = 1.7312480589;
			}
			else if (age == 48)
			{
				risk_5Y_average = 0.814639719;
				risk_10Y_average = 1.8166467243;
			}
			else if (age == 49)
			{
				risk_5Y_average = 0.8604567129;
				risk_10Y_average = 1.9002785226;
			}
			else if (age == 50)
			{
				risk_5Y_average = 0.9057802272;
				risk_10Y_average = 1.9817802881;
			}
			else if (age == 51)
			{
				risk_5Y_average = 0.9504300531;
				risk_10Y_average = 2.0607905074;
			}
			else if (age == 52)
			{
				risk_5Y_average = 0.9942261866;
				risk_10Y_average = 2.1369495822;
			}
			else if (age == 53)
			{
				risk_5Y_average = 1.0369888687;
				risk_10Y_average = 2.2099001348;
			}
			else if (age == 54)
			{
				risk_5Y_average = 1.0785386335;
				risk_10Y_average = 2.2792873637;
			}
			else if (age == 55)
			{
				risk_5Y_average = 1.1186963653;
				risk_10Y_average = 2.3447594549;
			}
			else if (age == 56)
			{
				risk_5Y_average = 1.1572833663;
				risk_10Y_average = 2.4059680578;
			}
			else if (age == 57)
			{
				risk_5Y_average = 1.1941214361;
				risk_10Y_average = 2.4625688323;
			}
			else if (age == 58)
			{
				risk_5Y_average = 1.2290329664;
				risk_10Y_average = 2.514222077;
			}
			else if (age == 59)
			{
				risk_5Y_average = 1.2618410501;
				risk_10Y_average = 2.5605934466;
			}
			else if (age == 60)
			{
				risk_5Y_average = 1.2923696103;
				risk_10Y_average = 2.6013547688;
			}
			else if (age == 61)
			{
				risk_5Y_average = 1.3204435495;
				risk_10Y_average = 2.63618497;
			}
			else if (age == 62)
			{
				risk_5Y_average = 1.345888922;
				risk_10Y_average = 2.6647711211;
			}
			else if (age == 63)
			{
				risk_5Y_average = 1.3685331334;
				risk_10Y_average = 2.6868096143;
			}
			else if (age == 64)
			{
				risk_5Y_average = 1.3882051698;
				risk_10Y_average = 2.7020074806;
			}
			else if (age == 65)
			{
				risk_5Y_average = 1.4047358598;
				risk_10Y_average = 2.7100838603;
			}
			else if (age == 66)
			{
				risk_5Y_average = 1.4179581739;
				risk_10Y_average = 2.7107716365;
			}
			else if (age == 67)
			{
				risk_5Y_average = 1.4277075651;
				risk_10Y_average = 2.7038192417;
			}
			else if (age == 68)
			{
				risk_5Y_average = 1.4338223551;
				risk_10Y_average = 2.6889926477;
			}
			else if (age == 69)
			{
				risk_5Y_average = 1.436144171;
				risk_10Y_average = 2.6660775451;
			}
			else if (age == 70)
			{
				risk_5Y_average = 1.4345184369;
				risk_10Y_average = 2.6348817206;
			}
			else if (age == 71)
			{
				risk_5Y_average = 1.4287949279;
				risk_10Y_average = 2.5952376349;
			}
			else if (age == 72)
			{
				risk_5Y_average = 1.4188283885;
				risk_10Y_average = 2.5470052021;
			}
			else if (age == 73)
			{
				risk_5Y_average = 1.4044792258;
				risk_10Y_average = 2.4900747681;
			}
			else if (age == 74)
			{
				risk_5Y_average = 1.3856142802;
				risk_10Y_average = 2.4243702792;
			}
		}
		else if (race == 5) // Hispanic
		{
			if (age == 35)
			{
				risk_5Y_average = 0.2487640896;
				risk_10Y_average = 0.6789848455;
			}
			else if (age == 36)
			{
				risk_5Y_average = 0.2829207566;
				risk_10Y_average = 0.7527102835;
			}
			else if (age == 37)
			{
				risk_5Y_average = 0.3185332961;
				risk_10Y_average = 0.8286164944;
			}
			else if (age == 38)
			{
				risk_5Y_average = 0.355472613;
				risk_10Y_average = 0.9064386875;
			}
			else if (age == 39)
			{
				risk_5Y_average = 0.3936095726;
				risk_10Y_average = 0.9859119724;
			}
			else if (age == 40)
			{
				risk_5Y_average = 0.4328150016;
				risk_10Y_average = 1.0667713451;
			}
			else if (age == 41)
			{
				risk_5Y_average = 0.4729596863;
				risk_10Y_average = 1.1487516637;
			}
			else if (age == 42)
			{
				risk_5Y_average = 0.5139143681;
				risk_10Y_average = 1.2315876147;
			}
			else if (age == 43)
			{
				risk_5Y_average = 0.5555497371;
				risk_10Y_average = 1.3150136693;
			}
			else if (age == 44)
			{
				risk_5Y_average = 0.5977364216;
				risk_10Y_average = 1.3987640328;
			}
			else if (age == 45)
			{
				risk_5Y_average = 0.6403449772;
				risk_10Y_average = 1.4825725847;
			}
			else if (age == 46)
			{
				risk_5Y_average = 0.6832458716;
				risk_10Y_average = 1.5661728148;
			}
			else if (age == 47)
			{
				risk_5Y_average = 0.7263094687;
				risk_10Y_average = 1.6492977534;
			}
			else if (age == 48)
			{
				risk_5Y_average = 0.7694060097;
				risk_10Y_average = 1.7316798992;
			}
			else if (age == 49)
			{
				risk_5Y_average = 0.8124055925;
				risk_10Y_average = 1.8130511469;
			}
			else if (age == 50)
			{
				risk_5Y_average = 0.8551781507;
				risk_10Y_average = 1.8931427165;
			}
			else if (age == 51)
			{
				risk_5Y_average = 0.8975934301;
				risk_10Y_average = 1.9716850885;
			}
			else if (age == 52)
			{
				risk_5Y_average = 0.9395209653;
				risk_10Y_average = 2.0484079481;
			}
			else if (age == 53)
			{
				risk_5Y_average = 0.9808300564;
				risk_10Y_average = 2.1230401419;
			}
			else if (age == 54)
			{
				risk_5Y_average = 1.0213897459;
				risk_10Y_average = 2.1953096546;
			}
			else if (age == 55)
			{
				risk_5Y_average = 1.061068797;
				risk_10Y_average = 2.264943609;
			}
			else if (age == 56)
			{
				risk_5Y_average = 1.0997356748;
				risk_10Y_average = 2.3316682985;
			}
			else if (age == 57)
			{
				risk_5Y_average = 1.1372585303;
				risk_10Y_average = 2.3952092576;
			}
			else if (age == 58)
			{
				risk_5Y_average = 1.1735051897;
				risk_10Y_average = 2.4552913828;
			}
			else if (age == 59)
			{
				risk_5Y_average = 1.2083431498;
				risk_10Y_average = 2.5116391112;
			}
			else if (age == 60)
			{
				risk_5Y_average = 1.2416395817;
				risk_10Y_average = 2.5639766714;
			}
			else if (age == 61)
			{
				risk_5Y_average = 1.2732613449;
				risk_10Y_average = 2.612028418;
			}
			else if (age == 62)
			{
				risk_5Y_average = 1.3030750142;
				risk_10Y_average = 2.6555192673;
			}
			else if (age == 63)
			{
				risk_5Y_average = 1.3309469217;
				risk_10Y_average = 2.69417525;
			}
			else if (age == 64)
			{
				risk_5Y_average = 1.3567432183;
				risk_10Y_average = 2.7277241997;
			}
			else if (age == 65)
			{
				risk_5Y_average = 1.3803299578;
				risk_10Y_average = 2.7558966001;
			}
			else if (age == 66)
			{
				risk_5Y_average = 1.4015732077;
				risk_10Y_average = 2.7784266124;
			}
			else if (age == 67)
			{
				risk_5Y_average = 1.4203391917;
				risk_10Y_average = 2.7950533091;
			}
			else if (age == 68)
			{
				risk_5Y_average = 1.4364944705;
				risk_10Y_average = 2.8055221412;
			}
			else if (age == 69)
			{
				risk_5Y_average = 1.4499061652;
				risk_10Y_average = 2.8095866669;
			}
			else if (age == 70)
			{
				risk_5Y_average = 1.4604422332;
				risk_10Y_average = 2.8070105735;
			}
			else if (age == 71)
			{
				risk_5Y_average = 1.4679718015;
				risk_10Y_average = 2.7975700206;
			}
			else if (age == 72)
			{
				risk_5Y_average = 1.4723655693;
				risk_10Y_average = 2.7810563379;
			}
			else if (age == 73)
			{
				risk_5Y_average = 1.4734962879;
				risk_10Y_average = 2.757279104;
			}
			else if (age == 74)
			{
				risk_5Y_average = 1.4712393301;
				risk_10Y_average = 2.7260696348;
			}
		}
		else	// use white for white, mixed, unknown
		{
			if (age == 35)
			{
				risk_5Y_average = 0.3148484523;
				risk_10Y_average = 0.8880099694;
			}
			else if (age == 36)
			{
				risk_5Y_average = 0.3622510142;
				risk_10Y_average = 0.9937820248;
			}
			else if (age == 37)
			{
				risk_5Y_average = 0.4125356074;
				risk_10Y_average = 1.1040124768;
			}
			else if (age == 38)
			{
				risk_5Y_average = 0.465480332;
				risk_10Y_average = 1.2182409184;
			}
			else if (age == 39)
			{
				risk_5Y_average = 0.5208631774;
				risk_10Y_average = 1.3360068991;
			}
			else if (age == 40)
			{
				risk_5Y_average = 0.5784620466;
				risk_10Y_average = 1.4568499938;
			}
			else if (age == 41)
			{
				risk_5Y_average = 0.6380547754;
				risk_10Y_average = 1.5803098473;
			}
			else if (age == 42)
			{
				risk_5Y_average = 0.6994191447;
				risk_10Y_average = 1.7059261978;
			}
			else if (age == 43)
			{
				risk_5Y_average = 0.7623328868;
				risk_10Y_average = 1.8332388784;
			}
			else if (age == 44)
			{
				risk_5Y_average = 0.8265736866;
				risk_10Y_average = 1.9617878004;
			}
			else if (age == 45)
			{
				risk_5Y_average = 0.8919191766;
				risk_10Y_average = 2.0911129209;
			}
			else if (age == 46)
			{
				risk_5Y_average = 0.9581469269;
				risk_10Y_average = 2.2207541966;
			}
			else if (age == 47)
			{
				risk_5Y_average = 1.0250344306;
				risk_10Y_average = 2.3502515283;
			}
			else if (age == 48)
			{
				risk_5Y_average = 1.0923590848;
				risk_10Y_average = 2.4791446999;
			}
			else if (age == 49)
			{
				risk_5Y_average = 1.1598981677;
				risk_10Y_average = 2.6069733168;
			}
			else if (age == 50)
			{
				risk_5Y_average = 1.2274288131;
				risk_10Y_average = 2.733276749;
			}
			else if (age == 51)
			{
				risk_5Y_average = 1.2947279824;
				risk_10Y_average = 2.857594086;
			}
			else if (age == 52)
			{
				risk_5Y_average = 1.3615724355;
				risk_10Y_average = 2.9794641108;
			}
			else if (age == 53)
			{
				risk_5Y_average = 1.4277387023;
				risk_10Y_average = 3.0984253014;
			}
			else if (age == 54)
			{
				risk_5Y_average = 1.4930030547;
				risk_10Y_average = 3.2140158704;
			}
			else if (age == 55)
			{
				risk_5Y_average = 1.5571414832;
				risk_10Y_average = 3.3257738539;
			}
			else if (age == 56)
			{
				risk_5Y_average = 1.6199296768;
				risk_10Y_average = 3.4332372617;
			}
			else if (age == 57)
			{
				risk_5Y_average = 1.6811430119;
				risk_10Y_average = 3.5359443055;
			}
			else if (age == 58)
			{
				risk_5Y_average = 1.7405565501;
				risk_10Y_average = 3.6334337199;
			}
			else if (age == 59)
			{
				risk_5Y_average = 1.7979450491;
				risk_10Y_average = 3.7252451953;
			}
			else if (age == 60)
			{
				risk_5Y_average = 1.8530829902;
				risk_10Y_average = 3.8109199438;
			}
			else if (age == 61)
			{
				risk_5Y_average = 1.905744626;
				risk_10Y_average = 3.8900014216;
			}
			else if (age == 62)
			{
				risk_5Y_average = 1.9557040525;
				risk_10Y_average = 3.9620362315;
			}
			else if (age == 63)
			{
				risk_5Y_average = 2.0027353124;
				risk_10Y_average = 4.026575236;
			}
			else if (age == 64)
			{
				risk_5Y_average = 2.0466125325;
				risk_10Y_average = 4.0831749095;
			}
			else if (age == 65)
			{
				risk_5Y_average = 2.0871101051;
				risk_10Y_average = 4.1313989628;
			}
			else if (age == 66)
			{
				risk_5Y_average = 2.1240029184;
				risk_10Y_average = 4.1708202751;
			}
			else if (age == 67)
			{
				risk_5Y_average = 2.1570666457;
				risk_10Y_average = 4.20102317;
			}
			else if (age == 68)
			{
				risk_5Y_average = 2.1860781025;
				risk_10Y_average = 4.2216060732;
			}
			else if (age == 69)
			{
				risk_5Y_average = 2.2108156811;
				risk_10Y_average = 4.23218459;
			}
			else if (age == 70)
			{
				risk_5Y_average = 2.2310598765;
				risk_10Y_average = 4.2323950416;
			}
			else if (age == 71)
			{
				risk_5Y_average = 2.2465939143;
				risk_10Y_average = 4.2218984948;
			}
			else if (age == 72)
			{
				risk_5Y_average = 2.257204496;
				risk_10Y_average = 4.2003853185;
			}
			else if (age == 73)
			{
				risk_5Y_average = 2.2626826773;
				risk_10Y_average = 4.1675802927;
			}
			else if (age == 74)
			{
				risk_5Y_average = 2.2628248957;
				risk_10Y_average = 4.123248287;
			}
		}
	}	
	return;
}

function create_result_page(age, race, relativeBC, biopsy, density)
{
	var race2 = race;
	if (race2 == 7)
	{
		race2 = 1;
	}
	calculate_risk (age, race2, relativeBC, biopsy, density);
	
	var header = document.getElementById("header");
	var headerText = header.firstChild;
	headerText.nodeValue = "Breast Cancer Risk Assessment Result";
	
	var tdForm = document.getElementById("tdRiskForm");	
	var riskForm = document.getElementById("riskForm");
	tdForm.removeChild(riskForm);
	
	var contentTable = document.getElementById("riskTable");	
	
	var resultSection = document.getElementById("divResult");
	resultSection.className = "result";
	
	var br1 = document.createElement("br");
	resultSection.appendChild(br1);
	
	// the patient's 5-year risk and 10-year risk
	var pRisk = document.createElement("p");
	pRisk.className = "result";
	resultSection.appendChild(pRisk);
	var textRisk = document.createTextNode("Based on the information provided, the woman's estimated risk for developing invasive breast cancer over the next 5 years is ");
	pRisk.appendChild(textRisk);
	var bRisk = document.createElement("b");
	pRisk.appendChild(bRisk);
	var textRisk2 = document.createTextNode(risk_5Y + "%");
	bRisk.appendChild(textRisk2);	
	var textRisk3 = document.createTextNode(", over the next 10 years is ");
	pRisk.appendChild(textRisk3);
	var bRisk2 = document.createElement("b");
	pRisk.appendChild(bRisk2);
	var textRisk4 = document.createTextNode(risk_10Y.toFixed(2) + "%");
	bRisk2.appendChild(textRisk4);
	var temp1 = document.createTextNode(".");
	pRisk.appendChild(temp1);
	
	// average risk for patients with same age group
	var displaySingleAverage = true;
	if (displaySingleAverage)
	{
		var pAvgRisk = document.createElement("p");
		pAvgRisk.className = "result";
		resultSection.appendChild(pAvgRisk);
		var textAvgRisk;
		// 5-year average
		if ((race == 1) || (race == 2) || (race == 3) || (race == 4) || (race == 5))
		{
			textAvgRisk = document.createTextNode("The average 5-year risk for a woman the same age and race/ethnicity is ");
		}
		else
		{
			textAvgRisk = document.createTextNode("The average risk 5-year for a white, non-Hispanic woman the same age is ");
		}
		pAvgRisk.appendChild(textAvgRisk);
		var bAvgRisk = document.createElement("b");
		pAvgRisk.appendChild(bAvgRisk);
		var textAvgRisk2 = document.createTextNode(risk_5Y_average.toFixed(2) + "%");
		bAvgRisk.appendChild(textAvgRisk2);
		var temp2 = document.createTextNode(". ");
		pAvgRisk.appendChild(temp2);
		// 10-year average
		if ((race == 1) || (race == 2) || (race == 3) || (race == 4) || (race == 5))
		{
			textAvgRisk = document.createTextNode("The average 10-year risk for a woman the same age and race/ethnicity is ");
		}
		else
		{
			textAvgRisk = document.createTextNode("The average 10-year risk for a white, non-Hispanic woman the same age is ");
		}
		pAvgRisk.appendChild(textAvgRisk);
		var bAvgRisk = document.createElement("b");
		pAvgRisk.appendChild(bAvgRisk);
		var textAvgRisk2 = document.createTextNode(risk_10Y_average.toFixed(2) + "%");
		bAvgRisk.appendChild(textAvgRisk2);
		var temp2 = document.createTextNode(".");
		pAvgRisk.appendChild(temp2);
	}
	else
	{
		// display average as a range with 2 numbers
		var pAvgRisk = document.createElement("p");
		pAvgRisk.className = "result";
		resultSection.appendChild(pAvgRisk);
		var textAvgRisk = document.createTextNode("The average risk for a woman within the same age group (" + get_age_group(age) + ") is between ");
		pAvgRisk.appendChild(textAvgRisk);
		var bAvgRiskLow = document.createElement("b");
		pAvgRisk.appendChild(bAvgRiskLow);
		var textAvgRiskLow = document.createTextNode(risk_avg_low.toFixed(2) + "%");
		bAvgRiskLow.appendChild(textAvgRiskLow);
		var and2 = document.createTextNode(" and ");
		pAvgRisk.appendChild(and2);
		var bAvgRiskHigh = document.createElement("b");
		pAvgRisk.appendChild(bAvgRiskHigh);
		var textAvgRiskHigh = document.createTextNode(risk_avg_high.toFixed(2) + "%");
		bAvgRiskHigh.appendChild(textAvgRiskHigh);
		var period2 = document.createTextNode(".");
		pAvgRisk.appendChild(period2);
	}	
		
	var pAnswers = document.createElement("p");
	pAnswers.className = "result";
	resultSection.appendChild(pAnswers);
	var textAnswer0 = document.createTextNode("These results are based upon the following answers about the woman:");
	pAnswers.appendChild(textAnswer0);
	var ulAnswers = document.createElement("ul");
	pAnswers.appendChild(ulAnswers);
	
	var liAnswer1 = document.createElement("li");
	liAnswer1.className = "result";
	ulAnswers.appendChild(liAnswer1);
	var textAnswer1 = document.createTextNode("Age: ");
	liAnswer1.appendChild(textAnswer1);
	var sAnswer1 = document.createElement("span");
	sAnswer1.className = "answer";
	liAnswer1.appendChild(sAnswer1);
	var textAnswer12 = document.createTextNode(age);
	sAnswer1.appendChild(textAnswer12);
	
	var liAnswer2 = document.createElement("li");
	liAnswer2.className = "result";
	ulAnswers.appendChild(liAnswer2);
	var textRace = document.createTextNode("Race/ethnicity: ");
	liAnswer2.appendChild(textRace);
	var sAnswer2 = document.createElement("span");
	sAnswer2.className = "answer";
	liAnswer2.appendChild(sAnswer2);
	var textAnswer2;
	if (race == 1)
	{
		textAnswer2 = document.createTextNode("White");
	}
	else if (race == 2)
	{
		textAnswer2 = document.createTextNode("Black");
	}
	else if (race == 3)
	{
		textAnswer2 = document.createTextNode("Asian");
	}
	else if (race == 4)
	{
		textAnswer2 = document.createTextNode("Native American");
	}
	else if (race == 5)
	{
		textAnswer2 = document.createTextNode("Hispanic");
	}
	else if (race == 6)
	{
		textAnswer2 = document.createTextNode("Other/Multiple races");
	}
	else // for cases when race is -1, 0 or 7
	{
		textAnswer2 = document.createTextNode("Unknown");
	}
	sAnswer2.appendChild(textAnswer2);
	
	if ((race < 1) || (race == 6) || (race == 7))
	{
		var textRace2 = document.createTextNode(" (data for white females is used to estimate the predicted risk)");
		liAnswer2.appendChild(textRace2);
	}
	
	var liAnswer3 = document.createElement("li");
	liAnswer3.className = "result";
	ulAnswers.appendChild(liAnswer3);
	var relativeText = document.createTextNode("First-degree relatives diagnosed with breast cancer: ");
	liAnswer3.appendChild(relativeText);
	var sAnswer3 = document.createElement("span");
	sAnswer3.className = "answer";
	liAnswer3.appendChild(sAnswer3);
	var textAnswer3;
	if (relativeBC == 0)
	{
		textAnswer3 = document.createTextNode("No");
	}
	if (relativeBC == 1)
	{
		textAnswer3 = document.createTextNode("Yes");
	}
	else if ((relativeBC > 1) || (relativeBC < 0))
	{
		textAnswer3 = document.createTextNode("Unknown");
	}
	sAnswer3.appendChild(textAnswer3);
	
	var liAnswer4 = document.createElement("li");
	liAnswer4.className = "result";
	ulAnswers.appendChild(liAnswer4);
	var biopsyText = document.createTextNode("Prior breast biopsy: ");
	liAnswer4.appendChild(biopsyText);
	var sAnswer4 = document.createElement("span");
	sAnswer4.className = "answer";
	liAnswer4.appendChild(sAnswer4);
	var textAnswer4;
	if (biopsy == 2)
	{
		textAnswer4 = document.createTextNode("None (no prior biopsy)");
	}
	else if (biopsy == 3)
	{
		textAnswer4 = document.createTextNode("Prior biopsy, unknown diagnosis");
	}
	else if (biopsy == 4)
	{
		textAnswer4 = document.createTextNode("Non-proliferative lesion");
	}
	else if (biopsy == 5)
	{
		textAnswer4 = document.createTextNode("Proliferative changes without atypia");
	}
	else if (biopsy == 6)
	{
		textAnswer4 = document.createTextNode("Proliferative changes with atypia");
	}
	else if (biopsy == 7)
	{
		textAnswer4 = document.createTextNode("Lobular carcinoma in situ");
	}
	else if ((biopsy > 7) || (biopsy < 2))
	{
		textAnswer4 = document.createTextNode("Unknown");
	}
	sAnswer4.appendChild(textAnswer4);
	
	var liAnswer5 = document.createElement("li");
	liAnswer5.className = "result";
	ulAnswers.appendChild(liAnswer5);
	var densityText = document.createTextNode("Breast density: ");
	liAnswer5.appendChild(densityText);
	var sAnswer5 = document.createElement("span");
	sAnswer5.className = "answer";
	liAnswer5.appendChild(sAnswer5);
	var textAnswer5;
	if (density == 1)
	{
		textAnswer5 = document.createTextNode("Almost entirely fatty");
	}
	else if (density == 2)
	{
		textAnswer5 = document.createTextNode("Scattered fibroglandular densities");
	}
	else if (density == 3)
	{
		textAnswer5 = document.createTextNode("Heterogeneously dense");
	}
	else
	{
		textAnswer5 = document.createTextNode("Extremely dense");
	}
	sAnswer5.appendChild(textAnswer5);
	
	var pReminder = document.createElement("p");
	pReminder.className = "remark";
	resultSection.appendChild(pReminder);
	var sReminder0 = document.createElement("span");
	sReminder0.className = "bold";
	pReminder.appendChild(sReminder0);
	var textReminder0 = document.createTextNode("Reminder:");
	sReminder0.appendChild(textReminder0);
	var textReminder1 = document.createTextNode(" The Breast Cancer Surveillance Consortium Risk Calculator was designed for use by health professionals. If you are not a health professional, you are encouraged to discuss the results and your personal risk of breast cancer with your health care provider.");
	pReminder.appendChild(textReminder1);
	
	var br2 = document.createElement("br");
	resultSection.appendChild(br2);
		
	var pAction = document.getElementById("pAction");
	var btnRestart = document.getElementById("mybutton");		
	pAction.removeChild(btnRestart);
	
	var aNewCalulation = document.createElement("a");
	aNewCalulation.className = "action";
	aNewCalulation.href = "calculator.htm";	
	var textNewCalculation = document.createTextNode("Start a New Risk Calculation");
	aNewCalulation.appendChild(textNewCalculation);	
	pAction.appendChild(aNewCalulation);
	
	window.scrollTo(0,0);
}



