package main

import (
	"flag"
	"fmt"
	"math/rand"
	"os"
	"strings"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/optimize"

	"github.com/kshedden/dstream/dstream"
	"github.com/kshedden/dstream/formula"
	"github.com/kshedden/duration"
	"github.com/kshedden/statmodel"
)

var (
	data dstream.Dstream

	nocenter = []string{"Time", "CvrgStart", "CvrgEnd", "HF", "HFDate", "DOB", "Weight"}
)

func drugGroupMain(vnames, ee []string) []string {
	for _, x := range vnames {
		if strings.HasPrefix(x, "TG") {
			// Omit these very rare categories
			f := false
			for _, y := range []int{5, 9, 12, 14, 18, 19, 22, 24} {
				if x == fmt.Sprintf("TG_%d", y) {
					f = true
				}
			}
			if !f {
				ee = append(ee, x)
			}
		}
	}
	return ee
}

func drugGroupInter(vnames, ee []string) []string {
	for _, x := range vnames {
		if strings.HasPrefix(x, "TG") {
			// Omit these very rare categories
			f := false
			for _, y := range []int{5, 9, 12, 14, 18, 19, 22, 24} {
				if x == fmt.Sprintf("TG_%d", y) {
					f = true
				}
			}
			if !f {
				ee = append(ee, x)
				ee = append(ee, "Age*"+x)
				ee = append(ee, "Female*"+x)
				ee = append(ee, "Age*Female*"+x)
			}
		}
	}
	return ee
}

func elixInter(vnames, ee []string) []string {
	for _, x := range vnames {
		if strings.HasPrefix(x, "Elix") && x != "Elix_CHF" {
			ee = append(ee, x)
			ee = append(ee, "Age*"+x)
			ee = append(ee, "Female*"+x)
			ee = append(ee, "Age*Female*"+x)
		}
	}
	return ee
}

func elixMain(vnames, ee []string) []string {
	for _, x := range vnames {
		if strings.HasPrefix(x, "Elix") && x != "Elix_CHF" {
			ee = append(ee, x)
		}
	}
	return ee
}

func procGroupMain(vnames, ee []string) []string {
	for _, x := range vnames {
		if strings.HasPrefix(x, "PG") {
			ee = append(ee, x)
		}
	}
	return ee
}

func procGroupInter(vnames, ee []string) []string {
	for _, x := range vnames {
		if strings.HasPrefix(x, "PG") {
			ee = append(ee, x)
			ee = append(ee, "Age*"+x)
			ee = append(ee, "Female*"+x)
			ee = append(ee, "Age*Female*"+x)
		}
	}
	return ee
}

func modeldata(model int) dstream.Dstream {

	var base = []string{"Age", "Female", "Age*Female"}
	var ee []string
	ee = append(ee, base...)
	vnames := data.Names()

	switch model {
	case 0:
		// Base model, demographics only
	case 1:
		ee = elixMain(vnames, ee)
	case 2:
		ee = elixInter(vnames, ee)
	case 3:
		ee = drugGroupMain(vnames, ee)
	case 4:
		ee = drugGroupInter(vnames, ee)
	case 5:
		ee = procGroupMain(vnames, ee)
	case 6:
		ee = procGroupInter(vnames, ee)
	case 7:
		ee = elixMain(vnames, ee)
		ee = drugGroupMain(vnames, ee)
		ee = procGroupMain(vnames, ee)
	case 8:
		ee = elixInter(vnames, ee)
		ee = drugGroupInter(vnames, ee)
		ee = procGroupInter(vnames, ee)
	}

	fml := strings.Join(ee, " + ")
	fmt.Printf(fml + "\n")

	keep := []string{"Time", "HF", "Weight"}
	dx := formula.New(fml, data).Keep(keep).Done()
	da := dstream.MemCopy(dx)

	return da
}

type rec struct {
	concordance float64
	l2w         float64
	result      *duration.PHResults
	kr          *statmodel.KnockoffResult
}

func ridge(mode int, l2w []float64, ko bool) error {

	da := modeldata(mode)

	if ko {
		// Names to knockoff
		var names []string
		for _, v := range da.Names() {
			if v != "Time" && v != "HF" && v != "Weight" {
				names = append(names, v)
			}
		}

		da.Reset()
		rand.Seed(323849)
		ko, err := statmodel.NewKnockoff(da, names)
		if err != nil {
			return err
		}
		ko.Reset()
		da = dstream.MemCopy(ko)
		fmt.Printf("lmin=%v\n", ko.CrossProdMinEig())
	}

	opt := optimize.DefaultSettings()
	opt.GradientThreshold = 1e-3

	rc := make(chan *rec, len(l2w))

	for _, w := range l2w {

		go func(w float64) {

			var l2wgt []float64
			for k := 0; k < len(da.Names())-3; k++ {
				l2wgt = append(l2wgt, w)
			}

			da.Reset()
			dx := dstream.Shallow(da)

			model := duration.NewPHReg(dx, "Time", "HF").OptSettings(opt).L2Weight(l2wgt).Weight("Weight")
			if !ko {
				// With knockoff we are already using normed data
				model = model.Norm()
			}
			model = model.Done()
			result, err := model.Fit()
			if err != nil {
				rc <- nil // Need to put something down the channel
				fmt.Printf("Error: %v\n", err)
				return
			}

			score := result.FittedValues(nil)

			dx.Reset()
			time := dstream.GetCol(dx, "Time").([]float64)
			dx.Reset()
			hf := dstream.GetCol(dx, "HF").([]float64)
			c := duration.NewConcordance(time, hf, score).Done()
			var kr *statmodel.KnockoffResult
			if ko {
				kr = statmodel.NewKnockoffResult(result, false)
			}
			rc <- &rec{c.Concordance(365), w, result, kr}
		}(w)
	}

	fname := fmt.Sprintf("coeff_%d.txt", mode)
	if ko {
		fname = strings.Replace(fname, ".txt", "_ko.txt", 1)
	}
	fid, err := os.Create(fname)
	if err != nil {
		panic(err)
	}
	defer fid.Close()
	for k := 0; k < len(l2w); k++ {
		r := <-rc
		if r != nil {
			fid.Write([]byte(fmt.Sprintf("L2=%f  concordance=%f\n", r.l2w, r.concordance)))
			if r.kr != nil {
				names := r.kr.Names()
				params := r.kr.Params()
				stat := r.kr.Stat()
				fdr := r.kr.FDR()
				for i := range names {
					fid.WriteString(fmt.Sprintf("%s\t%f\t%f\t%f\n", names[i], params[i], stat[i], fdr[i]))
				}
			} else {
				fid.WriteString(r.result.Summary() + "\n")
			}
		}
	}

	return nil
}

func center(data dstream.Dstream) dstream.Dstream {

	means := make([]float64, data.NumVar())
	ntot := 0

	nc := make(map[string]bool)
	for _, na := range nocenter {
		nc[na] = true
	}

	data.Reset()
	for data.Next() {
		for k := range data.Names() {
			z := data.GetPos(k).([]float64)
			if k == 0 {
				ntot += len(z)
			}
			means[k] += floats.Sum(z)
		}
	}
	for k := range means {
		means[k] /= float64(ntot)
	}

	for k, na := range data.Names() {
		if nc[na] {
			continue
		}

		mn := means[k]
		f := func(x interface{}) {
			z := x.([]float64)
			for i := range z {
				z[i] -= mn
			}
		}

		data = dstream.Mutate(data, na, f)
	}

	return data
}

func genvars() {

	f := func(v map[string]interface{}, x interface{}) {
		wt := x.([]float64)
		hf := v["HF"].([]float64)
		for i := range hf {
			if hf[i] == 1 {
				wt[i] = 1
			} else {
				wt[i] = 10
			}
		}
	}
	data = dstream.Generate(data, "Weight", f, "float64")

	f = func(v map[string]interface{}, x interface{}) {
		age := x.([]float64)
		dob := v["DOB"].([]float64)
		cvrgstart := v["CvrgStart"].([]float64)
		for i := range dob {
			y := 1960 + (cvrgstart[i]-1960)/365.25
			age[i] = y - dob[i]
		}
	}
	data = dstream.Generate(data, "Age", f, "float64")

	data.Reset()
	st := dstream.Describe(data)
	for k, v := range st {
		fmt.Printf("%s\n", k)
		fmt.Printf("%+v\n", v)
	}
}

func main() {

	var ko bool
	flag.BoolVar(&ko, "knockoff", false, "Use knockoff method")
	flag.Parse()

	data = dstream.NewBCols("data", 100000).Done()
	genvars()
	data = center(data)

	stats := dstream.Describe(data)
	for k, v := range stats {
		fmt.Printf("%-20s %f\n", k, v.SD)
	}

	l2w := []float64{0.05, 0.1, 0.2, 0.4, 0.8, 1.6}
	for k := 0; k < 9; k++ {
		err := ridge(k, l2w, ko)
		if err != nil {
			print(err)
		}
	}
}
