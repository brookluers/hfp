package utils

// Drec describes one subject in the data set.
type Drec struct {

	// Indicator that the subject has heart failure
	Hf bool

	// Date at which the subject first had heart failure
	HfDate uint16

	// First date of coverage
	CvrgStart uint16

	// Last date of coverage
	CvrgEnd uint16

	// Date of birth
	DOB uint16

	// Sex
	Sex uint8

	// Array of Elixhauser indicators
	Elix []int

	// Array of drug therapeutic group indicators
	Thrgrp []int

	// Array of procedure group codes
	Procgrp []int
}
