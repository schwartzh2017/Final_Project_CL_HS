// Copyright 2012-2015, Breast Cancer Surveillance Consortium
// Developed by Jambeyang Research, LLC., http://www.jambeyang.com, email: wlin@jambeyang.com
// Last modified on 4/27/2015

function checkPriorBC() {
	if (document.riskForm.priorBC.selectedIndex == 1)
	{
		alert("This calculation is not applicable to women with any previous diagnosis of breast cancer, DCIS, prior breast augmentation, or prior mastectomy.");
	}
	else if (document.riskForm.priorBC.selectedIndex < 1)
	{
		alert("Please answer the question about the woman's prior breast diseases (Question 1).");
	}
	else
	{
		return true;
	}
	return false;
}
		
function checkAge()
{
	if (document.riskForm.age.selectedIndex == 1)
	{
		alert("This calculation is not applicable to women younger than age 35.");
	}
	else if (document.riskForm.age.selectedIndex == (document.riskForm.age.length - 1))
	{
		alert("This calculation is not applicable to women older than age 74.");
	}
	else if (document.riskForm.age.selectedIndex < 1)
	{
		alert("Please provide the woman's age (Question 2).");
	}
	else
	{
		return true;
	}
	return false;
}

function checkRace()
{
	if (document.riskForm.race.selectedIndex < 1)
	{
		alert("Please provide the woman's race/ethnicity (Question 3).");
		return false;
	}
	else
	{
		return true;
	}
}
		
function checkRelativeBC()
{
	if (document.riskForm.relativeBC.selectedIndex < 1)
	{
		alert("Please answer the question about the woman's first-degree relatives (Question 4).");
		return false;
	}
	else
	{
		return true;
	}
}
		
function checkBiopsy()
{
	if (document.riskForm.biopsy.selectedIndex < 1)
	{
		alert("Please answer the question about the woman's prior breast biopsies (Question 5).");
		return false;
	}
	else
	{
		return true;
	}
}
		
function checkDensity()
{
	if (document.riskForm.density.selectedIndex < 1)
	{
		alert("Please provide the woman's breast density (Question 6).");
		return false;
	}
	else
	{
		return true;
	}
}