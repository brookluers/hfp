/*
Create binary columns from the variables in hfdat.gob.gz.
*/

package main

import (
	"compress/gzip"
	"encoding/binary"
	"encoding/gob"
	"encoding/json"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"path"
	"sort"
	"strings"

	"github.com/brookluers/hfp/utils"
)

// xws is a structure to manage output processing for one variable.
type xws struct {

	// Write directly to the file
	fw io.WriteCloser

	// Write compressed data to the file
	zw io.WriteCloser

	xt func(*utils.Drec) float64
}

// Init sets up an xws value to write t the given file, using values extracted from the
// data record using the given extractor function.
func (xw *xws) Init(fname string, f func(*utils.Drec) float64) {

	var err error
	xw.fw, err = os.Create(path.Join("data", fmt.Sprintf("%s.bin.gz", fname)))
	if err != nil {
		panic(err)
	}

	xw.zw = gzip.NewWriter(xw.fw)
	xw.xt = f
}

// Add extracts and writes one value from the record to the binary column file.
func (xw *xws) Add(r *utils.Drec) {
	err := binary.Write(xw.zw, binary.LittleEndian, xw.xt(r))
	if err != nil {
		panic(err)
	}
}

// Close closes the io writers.
func (xw *xws) Close() {
	xw.zw.Close() // order is important here
	xw.fw.Close()
}

func maindata() {

	// Set up for reading the gob.
	f, err := os.Open("hfdat.gob.gz")
	if err != nil {
		panic(err)
	}
	defer f.Close()
	g, err := gzip.NewReader(f)
	if err != nil {
		panic(err)
	}
	defer g.Close()

	dec := gob.NewDecoder(g)

	// Elixhauser category names
	var elxn []string
	err = dec.Decode(&elxn)
	if err != nil {
		panic(err)
	}

	// Set up for writing the Elixhauser values.
	elxw := make([]*xws, len(elxn))
	for i := range elxw {
		el := new(xws)
		name := fmt.Sprintf("Elix_%s", elxn[i])
		ii := i
		el.Init(name,
			func(r *utils.Drec) float64 {
				j := sort.SearchInts(r.Elix, ii)
				if j == len(r.Elix) || r.Elix[j] != ii {
					return 0
				}
				return 1
			},
		)
		elxw[i] = el
		defer el.Close()
	}

	// Set up for writing the drug therapeutic group values
	tgw := make([]*xws, 31)
	for i := 0; i < 31; i++ {
		tg := new(xws)
		name := fmt.Sprintf("TG_%02d", i)
		ii := i
		tg.Init(name,
			func(r *utils.Drec) float64 {
				j := sort.SearchInts(r.Thrgrp, ii)
				if j == len(r.Thrgrp) || r.Thrgrp[j] != ii {
					return 0
				}
				return 1
			},
		)
		tgw[i] = tg
		defer tg.Close()
	}

	hfdate := new(xws)
	hfdate.Init("HFDate", func(r *utils.Drec) float64 { return float64(r.HfDate) })
	defer hfdate.Close()

	hf := new(xws)
	hf.Init(
		"HF",
		func(r *utils.Drec) float64 {
			if r.Hf {
				return 1
			}
			return 0
		},
	)
	defer hf.Close()

	dob := new(xws)
	dob.Init("DOB", func(r *utils.Drec) float64 { return float64(r.DOB) })
	defer dob.Close()

	cvrgstart := new(xws)
	cvrgstart.Init("CvrgStart", func(r *utils.Drec) float64 { return float64(r.CvrgStart) })
	defer cvrgstart.Close()

	cvrgend := new(xws)
	cvrgend.Init("CvrgEnd", func(r *utils.Drec) float64 { return float64(r.CvrgEnd) })
	defer cvrgend.Close()

	female := new(xws)
	female.Init("Female", func(r *utils.Drec) float64 { return float64(r.Sex - 1) })
	defer female.Close()

	tim := new(xws)
	tim.Init("Time",
		func(r *utils.Drec) float64 {
			if r.Hf {
				return float64(r.HfDate - r.CvrgStart - 365)
			}
			return float64(r.CvrgEnd - r.CvrgStart - 365)
		},
	)
	defer tim.Close()

	nrec := 0
	for {
		var r utils.Drec
		err := dec.Decode(&r)
		if err == io.EOF {
			break
		} else if err != nil {
			panic(err)
		}

		// Needs to match selection in reduce.go.
		if r.DOB > 1970 {
			continue
		}
		nrec++

		hfdate.Add(&r)
		hf.Add(&r)
		dob.Add(&r)
		cvrgstart.Add(&r)
		cvrgend.Add(&r)
		tim.Add(&r)
		female.Add(&r)

		for _, e := range elxw {
			e.Add(&r)
		}

		for _, t := range tgw {
			t.Add(&r)
		}
	}

	fmt.Printf("Processed %d records\n", nrec)
}

// dtypes updates the dtype information in the data directory.
func dtypes() {

	dt := make(map[string]string)

	fi, err := ioutil.ReadDir("data")
	if err != nil {
		panic(err)
	}

	for _, f := range fi {
		a := f.Name()
		if strings.HasSuffix(a, ".bin.gz") {
			a = strings.Replace(a, ".bin.gz", "", -1)
			dt[a] = "float64"
		}
	}

	out, err := os.Create("data/dtypes.json")
	if err != nil {
		panic(err)
	}
	enc := json.NewEncoder(out)
	err = enc.Encode(&dt)
	if err != nil {
		panic(err)
	}
}

func main() {

	maindata()
	dtypes()
}
