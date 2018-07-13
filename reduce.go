/*
Reduce a large spares matrix to factors, using an approximate SVD.

Right now this is hard-coded to work on the Procgrp data from hfdat.gob.gz, but
could be generalized to work on other data.

The resulting factors are stored in the 'data' directory as binary columns.
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
	"math"
	"os"
	"path"
	"strings"

	"gonum.org/v1/gonum/mat"

	"github.com/brookluers/dimred"
	"github.com/brookluers/hfp/utils"
)

func doFactorize(nfac, npow int) (int, *mat.Dense, *mat.Dense) {

	// Setup gob file reader
	inf, err := os.Open("hfdat.gob.gz")
	if err != nil {
		panic(err)
	}
	defer inf.Close()
	gnf, err := gzip.NewReader(inf)
	if err != nil {
		panic(err)
	}
	dec := gob.NewDecoder(gnf)

	// The first element is the Elixhauser category names
	var elxcat []string
	err = dec.Decode(&elxcat)
	if err != nil {
		panic(err)
	}

	// The sparse matrix is represented as mat[row[i], col[i]] = dat[i]
	var row, col []int
	var dat []float64

	var nrec int
	for {

		// Read one subject record
		var r utils.Drec
		err := dec.Decode(&r)
		if err == io.EOF {
			break
		} else if err != nil {
			panic(err)
		}

		// We will eventually select ages 50-65 so this is
		// just a pre-filter.  This selection needs to match
		// the analogous selection in data.go.
		if r.DOB > 1970 {
			continue
		}
		nrec++

		// Add an entry to the sparse matrix
		for _, g := range r.Procgrp {
			row = append(row, nrec)
			col = append(col, g)
			dat = append(dat, 1)
		}
	}
	fmt.Printf("Processsed %d records\n", nrec)

	// Run the approximate SVD
	spm := dimred.NewSPM(row, col, dat, nrec, 500)
	sv := new(dimred.RSVD)
	sv.Factorize(spm, nfac, npow)
	umat := sv.UTo(nil)
	vmat := sv.VTo(nil)

	fmt.Printf("values: %f\n", sv.Values(nil))

	// Center and scale the columns
	nrow, ncol := umat.Dims()
	for j := 0; j < ncol; j++ {
		v := umat.ColView(j)

		mn := 0.0
		for i := 0; i < nrow; i++ {
			mn += v.AtVec(i)
		}
		mn /= float64(nrow)

		for i := 0; i < nrow; i++ {
			umat.Set(i, j, umat.At(i, j)-mn)
		}

		sc := 0.0
		for i := 0; i < nrow; i++ {
			z := v.AtVec(i)
			sc += z * z
		}
		sc = math.Sqrt(sc)

		for i := 0; i < nrow; i++ {
			umat.Set(i, j, umat.At(i, j)/sc)
		}
	}

	return nrec, umat, vmat
}

//  store writes column in binary form
func store(ma *mat.Dense, dir, pre string, sf float64) {

	var out []io.WriteCloser
	nrow, ncol := ma.Dims()

	// Create the destination file writers
	for j := 0; j < ncol; j++ {
		fname := path.Join(dir, fmt.Sprintf("%s_%03d.bin.gz", pre, j))
		f, err := os.Create(fname)
		if err != nil {
			panic(err)
		}
		defer f.Close()
		z := gzip.NewWriter(f)
		defer z.Close()
		out = append(out, z)
	}

	// Write the data
	for i := 0; i < nrow; i++ {
		for j := 0; j < ncol; j++ {
			err := binary.Write(out[j], binary.LittleEndian, sf*ma.At(i, j))
			if err != nil {
				panic(err)
			}
		}
	}
}

// storev stores the v matrix of the approximate SVD in a file.
func storev(m *mat.Dense, fname string) {

	f, err := os.Create(fname)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	g := gzip.NewWriter(f)
	defer g.Close()

	m.MarshalBinaryTo(g)
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

	// Number of factors to extract
	nfac := 20

	// Number of power iterations to apply during the approximate SVD
	npow := 5

	n, umat, vmat := doFactorize(nfac, npow)

	store(umat, "data", "PG", math.Sqrt(float64(n)))

	storev(vmat, "procgrp_v.bin.gz")

	dtypes()
}
