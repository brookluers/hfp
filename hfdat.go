/*
Create a gob file containing records for all eligible people with heart failure.
Also include a random sample of eligible people without heart failure.
*/

package main

import (
	"compress/gzip"
	"encoding/binary"
	"encoding/gob"
	"encoding/json"
	"fmt"
	"hash/crc32"
	"io"
	"log"
	"os"
	"sort"
	"strings"

	"github.com/kshedden/dstream/dstream"
	"github.com/kshedden/gocols/config"
	"github.com/kshedden/hfp/utils"
)

const (
	// Chunk size for reading raw data
	csize = 100000

	// Exclude a person for one year if they do not have this many
	// days of coverage in the year
	mincoverage = 360

	// Process this number of buckets in parallel
	concurrency = 100
)

var (
	// Source file names (A, I, etc.)
	dix []string

	// Map from ICD codes to our integer representation
	dxcodes map[string]int

	// Int codes corresponding to heart faiure ICD codes
	hfcodes []int

	// Names of Elixhauser categories
	elxcat []string

	// Int codes corresponding to ICD codes for each Elixhauser category
	elix [][]int

	logger *log.Logger

	oif io.WriteCloser
	oig io.WriteCloser

	rslt chan utils.Drec
	sem  chan bool

	// Configuration files
	aconf, oconf, sconf, iconf, fconf, dconf *config.Config

	// Driectory paths to the configuration files
	adir, odir, sdir, idir, fdir, ddir string
)

// bcomp compresses a Boolean array into a list of indices where the
// values are true.  For example, [true, false, true] becomes [0, 2].
func bcomp(x []bool) []int {
	var y []int
	for i, v := range x {
		if v {
			y = append(y, i)
		}
	}
	return y
}

// convertICDCodes takes a list of ICD9 codes in string form and maps
// them to their corresponding integer codes.
func convertICDCodes(r []string) []int {

	var i []int
	for _, x := range r {
		c, ok := dxcodes[x]
		if ok {
			i = append(i, c)
		}
	}

	return i
}

// matchset returns true or false based on whether y is an element of
// the sorted array x.
func matchset(x []int, y int) bool {

	j := sort.SearchInts(x, y)

	if j == len(x) || x[j] != y {
		return false
	}

	return true
}

// readRawElix loads and returns the Elixhauser ICD codes for either
// version 9 or version 10.
func readRawElix(ver int) map[string][]string {

	fn := fmt.Sprintf("elix%d.json", ver)
	fid, err := os.Open(fn)
	if err != nil {
		panic(err)
	}
	defer fid.Close()

	dec := json.NewDecoder(fid)
	elx := make(map[string][]string)
	err = dec.Decode(&elx)
	if err != nil {
		panic(err)
	}

	return elx
}

// getElix returns an array of arrays containing sorted, distinct codes for each Elixhauser category,
// and a list of names for the categories.
func getElix() ([][]int, []string) {

	var elxcat []string

	// The ICD9 and ICD10 codes for Elixhauser categories
	elix9 := readRawElix(9)
	elix10 := readRawElix(10)

	// Get the Elixhauser category names (assume Elix9 and Elix10 category
	// names are the same).
	for k := range elix9 {
		elxcat = append(elxcat, k)
	}
	sort.Sort(sort.StringSlice(elxcat))

	// Get the integer codes for all Dx codes in each Elixhauser
	// category, mixing ICD9 and ICD10 codes.
	var elx [][]int
	for _, cat := range elxcat {

		var u []int
		u = append(u, convertICDCodes(elix9[cat])...)
		u = append(u, convertICDCodes(elix10[cat])...)

		// Sort and remove duplicates
		sort.Sort(sort.IntSlice(u))
		j := 0
		for i := range u {
			if i == 0 || u[i] != u[i-1] {
				u[j] = u[i]
				j++
			}
		}
		u = u[0:j]
		elx = append(elx, u)
	}

	return elx, elxcat
}

// edays takes a four digit year, and returns the elapsed days from 1-1-1960.
func edays(year int) uint16 {
	return uint16(365.25 * (float64(year) - 1960))
}

// setup Bucket returns an array of dstreams corresponding to
// A, O, S, I, F and D tables.  Relevant columns are selected
// from each table.
func setupBucket(bk int) []dstream.Dstream {

	var data []dstream.Dstream

	for _, v := range []struct {
		dirpath string
		include []string
	}{
		{
			// A tables
			dirpath: adir,
			include: []string{"Enrolid", "Year", "Memdays", "Dobyr", "Region", "Emprel", "Sex"},
		},
		{
			// O tables
			dirpath: odir,
			include: []string{"Enrolid", "Svcdate", "Dx1", "Dx2", "Dx3", "Dx4"},
		},
		{
			// S tables
			dirpath: sdir,
			include: []string{"Enrolid", "Svcdate", "Dx1", "Dx2"},
		},
		{
			// I tables
			dirpath: idir,
			include: []string{"Enrolid", "Admdate", "Dx1", "Dx2", "Dx3", "Dx4", "Dx5",
				"Dx6", "Dx7", "Dx8", "Dx9", "Dx10", "Dx11", "Dx12", "Dx13",
				"Dx14", "Dx15"},
		},
		{
			// F tables
			dirpath: fdir,
			include: []string{"Enrolid", "Svcdate", "Dx1", "Dx2", "Dx3", "Dx4", "Dx5",
				"Dx6", "Dx7", "Dx8", "Dx9"},
		},
		{
			// D tables
			dirpath: ddir,
			include: []string{"Enrolid", "Svcdate", "Thergrp"},
		},
	} {
		bp := config.BucketPath(bk, v.dirpath)
		da := dstream.NewBCols(bp, csize).Include(v.include).Done()
		da = dstream.Segment(da, []string{"Enrolid"})
		data = append(data, da)
	}

	return data
}

// getEligible finds the longest consecutive series of years during which the subject is eligible.
func getEligible(ayear []uint16, amemdays []uint16) (int, int) {

	// Years with full coverage
	ym := make(map[int]bool)
	for k := range ayear {
		if amemdays[k] >= mincoverage {
			ym[int(ayear[k])] = true
		}
	}

	if len(ym) == 0 {
		return -1, -1
	}

	year0, year1 := -1, -1
	mxd := 0
	for y := range ym {
		yy := y + 1
		for ym[yy] {
			// advance until the year is missing
			yy++
		}
		if yy-y > mxd {
			mxd = yy - y
			year0 = y
			year1 = yy
		}
	}

	return year0, year1
}

// dxPos returns the positions of all Dx variables in each of the tables.
func dxPos(data []dstream.Dstream) [][]int {

	var dpx [][]int
	for j := range dix {
		vpos := dstream.VarPos(data[j])
		var u []int
		for k, v := range vpos {
			if strings.HasPrefix(k, "Dx") {
				u = append(u, v)
			}
		}
		dpx = append(dpx, u)
	}

	return dpx
}

func dobucket(k int) {

	defer func() { <-sem }()

	logger.Printf("Starting bucket %d\n", k)

	ha := crc32.NewIEEE()

	data := setupBucket(k)
	adata := data[0]

	dpx := dxPos(data)

	// Positions of the date variable in each file
	datename := []string{"Year", "Svcdate", "Svcdate", "Admdate", "Svcdate", "Svcdate"}

	wk := dstream.NewJoin(data, []string{"Enrolid", "Enrolid", "Enrolid", "Enrolid",
		"Enrolid", "Enrolid"})

	// Loop over subjects
	for js := 0; wk.Next(); js++ {

		if js%100000 == 0 {
			logger.Printf("Bucket %d: %d", k, js)
		}

		// Get the eligible years.
		ayear := adata.Get("Year").([]uint16)
		amemdays := adata.Get("Memdays").([]uint16)
		year0, year1 := getEligible(ayear, amemdays)
		if year0 == -1 || year1-year0 == 1 {
			// Note enough coverage
			continue
		}

		// First/last day of coverage window
		d0 := edays(year0)
		d1 := edays(year1)

		aenrolid := adata.Get("Enrolid").([]uint64)
		adobyr := adata.Get("Dobyr").([]uint16)
		asex := adata.Get("Sex").([]uint8)

		// Heart failure status
		var hf = false
		var hfdate uint16

		elx := make([]bool, len(elix)) // Elixhauser categories
		thg := make([]bool, 31)        // Therapeutic groups for drugs
		pcx := make([]bool, 500)       // Procedure grops

		// Loop over tables
		for j, vb := range dix {

			// Skip files that contain no useful information
			if vb == "A" || !wk.Status[j] {
				continue
			}

			// The time values
			sd := data[j].Get(datename[j]).([]uint16)

			// Drug information
			if vb == "D" {
				thr := data[j].Get("Thergrp").([]uint8)
				for i, t := range thr {
					if sd[i] >= d0 && sd[i] <= d0+365 && t > 0 && t <= 31 {
						thg[t-1] = true
					}
				}
			}

			// Procgrp information
			if vb == "O" {
				pcg := data[j].Get("Procgrp").([]uint16)
				for i, p := range pcg {
					if sd[i] >= d0 && sd[i] <= d0+365 && p > 0 {
						pcx[p-1] = true
					}
				}
			}

			// Loop over the Dx variables in each file
			for _, k := range dpx[j] {

				// Loop through time
				for i, y := range data[j].GetPos(k).([]uint64) {

					// Check for heart failure
					if matchset(hfcodes, int(y)) {
						if !hf {
							// First observation of an HF code for this subject
							hfdate = sd[i]
						} else {
							if sd[i] < hfdate {
								// Earliest HF code for this subject
								hfdate = sd[i]
							}
						}
						hf = true
					}

					// Update Elixhauser for first year after entry.
					if sd[i] >= d0 && sd[i] <= d0+365 {
						for q, ex := range elix {
							elx[q] = elx[q] || matchset(ex, int(y))
						}
					}
				}
			}
		}

		// Heart failure within first year, exclude
		if hf && (hfdate < d0+365) {
			continue
		}

		// Take a random subsample of non-cases
		if !hf {
			binary.Write(ha, binary.LittleEndian, aenrolid[0])
			if ha.Sum32()%10 != 0 {
				continue
			}
		}

		// A record for this subject
		r := utils.Drec{
			Hf:        hf,
			HfDate:    hfdate,
			CvrgStart: d0,
			CvrgEnd:   d1,
			Elix:      bcomp(elx),
			Thrgrp:    bcomp(thg),
			Procgrp:   bcomp(pcx),
			DOB:       adobyr[0],
			Sex:       asex[0],
		}

		rslt <- r
	}
}

func setuplog() {
	fid, err := os.Create("hfdat.log")
	if err != nil {
		panic(err)
	}
	logger = log.New(fid, "", log.Ltime)
}

// harvest collects data from all buckets
func harvest() {

	enc := gob.NewEncoder(oig)

	// First write the Elixhauser category names to the gob.
	err := enc.Encode(elxcat)
	if err != nil {
		panic(err)
	}

	// The remainder of the gob file is the sequence of records, 1 per subject.
	for r := range rslt {
		err := enc.Encode(r)
		if err != nil {
			panic(err)
		}
	}
}

// setHF finds the HF ICD codes.
func setHF() {

	ii := -1
	for j, c := range elxcat {
		if c == "CHF" {
			ii = j
			break
		}
	}
	if ii == -1 {
		panic("can't find CHF category\n")
	}
	hfcodes = elix[ii]
}

func main() {

	if len(os.Args) != 7 {
		os.Stderr.WriteString("Usage:\nhfdat aconfig oconfig sconfig iconfig fconfig dconfig\n")
		os.Exit(1)
	}

	setuplog()

	var err error
	oif, err = os.Create("hfdat.gob.gz")
	if err != nil {
		panic(err)
	}
	defer oif.Close()

	oig = gzip.NewWriter(oif)
	defer oig.Close()

	adir = os.Args[1]
	odir = os.Args[2]
	sdir = os.Args[3]
	idir = os.Args[4]
	fdir = os.Args[5]
	ddir = os.Args[6]

	aconf = config.GetConfig(adir)
	oconf = config.GetConfig(odir)
	sconf = config.GetConfig(sdir)
	iconf = config.GetConfig(idir)
	fconf = config.GetConfig(fdir)
	dconf = config.GetConfig(ddir)

	dxcodes = config.GetFactorCodes("Dx", oconf)

	dix = []string{"A", "O", "S", "I", "F", "D"}

	sem = make(chan bool, concurrency)
	rslt = make(chan utils.Drec, 200)

	elix, elxcat = getElix()

	setHF()

	go harvest()

	for k := 0; k < int(oconf.NumBuckets); k++ {
		sem <- true
		go dobucket(k)
	}

	for k := 0; k < concurrency; k++ {
		sem <- true
	}

	close(rslt)
}
