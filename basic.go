package main

import (
	"flag"
	"fmt"
	"math/rand"
	"os"
	"strings"
	"time"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/optimize"

	"github.com/brookluers/dimred"
	"github.com/brookluers/dstream/dstream"
	"github.com/brookluers/dstream/formula"
	"github.com/brookluers/duration"
	"github.com/brookluers/statmodel"
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
			// BGL (6/21/18) TG_11 and TG_23 had zero variance, try excluding
			// for _, y := range []int{5, 9, 11, 12, 14, 18, 19, 22, 23, 24} {
			//	if x == fmt.Sprintf("TG_%d", y) {
			//		f = true
			//	}
			//}
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
			// BGL (6/21/18) TG_11 and TG_23 had zero variance, try excluding
			//for _, y := range []int{5, 9, 11, 12, 14, 18, 19, 22, 23, 24} {
			//	if x == fmt.Sprintf("TG_%d", y) {
			//		f = true
			//	}
			//}
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

func modeldata(model int, fl_save bool, fl_fullrank bool, fl_qr bool, ko bool) dstream.Dstream {

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
	// keep variables in dstream but not in formula
	dx := formula.New(fml, data).Keep(keep).Done()

	fmt.Printf("\n---Variable names after formula parsing---%v\n\n", dx.Names())

	if fl_fullrank {
		lfn := "log_fr.txt"
		if fl_qr {
			lfn = strings.Replace(lfn, ".txt", "_qr.txt", 1)
		}

		if ko {
			lfn = strings.Replace(lfn, ".txt", "_ko.txt", 1)
		}
		fmt.Printf("log file name = %v\n", lfn)
		var pcheck int
		var frcpr []float64
		if fl_qr {
			fmt.Printf("\n--Finding set of linearly independent columns using rank-revealing QR--\n")
			fr := dimred.NewRRQR(dx).Keep("Time", "HF", "Weight").LogFile(lfn).Tol(1e-5).Done()
			pcheck = fr.DimCheck()
			frcpr = fr.CPR()
			dx = fr.Data()
		} else {
			fmt.Printf("\n--Finding set of linearly independent columns using Cholesky--\n")
			fr := dimred.NewFullRank(dx).Keep("Time", "HF", "Weight").LogFile(lfn).Tol(1e-5).Done()
			pcheck = fr.DimCheck()
			frcpr = fr.CPR()
			dx = fr.Data()
		}

		gmat := mat.NewDense(pcheck, pcheck, frcpr)
		fmt.Printf("Upper corner of Gram matrix: %v\n\n",
			mat.Formatted(gmat, mat.Prefix(" "), mat.Excerpt(3)))

		var gsvd mat.SVD
		ok := gsvd.Factorize(gmat, mat.SVDNone)
		if !ok {
			fmt.Printf("Failed to compute SVD of gram matrix\n")
		} else {
			fmt.Printf("\nSingular values of gram matrix: %v\n", gsvd.Values(nil))
		}

		fmt.Printf("\nNames of full-rank column set: \n%v\n", dx.Names())
		fmt.Printf("\nNumber of (full-rank) columns: %v\n", len(dx.Names()))
	}

	da := dstream.MemCopy(dx)

	if fl_save {
		nrec := da.NumObs()
		if nrec > 0 {
			rand.Seed(718191)
			rsavefunc := func(v map[string]interface{}, x interface{}) {
				saveit := x.([]float64)
				hf := v["HF"].([]float64)
				for i := range hf {
				     saveit[i] = rand.Float64()
				}
			}
			saveProp := 0.05
			filterSave := func(x interface{}, b []bool) bool {
    				    v := x.([]float64)
    				      var any bool
				      for i := range v {
			                  b[i] = v[i] < saveProp
        			          any = any || !b[i]
    				      }
				      return any
		        }
			
			da2 := dstream.Generate(da, "SaveRand", rsavefunc, "float64")
			da2 = dstream.Filter(da2, map[string]dstream.FilterFunc{"SaveRand": filterSave})
			
			werr := dstream.ToCSV(da2).Filename(fmt.Sprintf("/scratch/stats_flux/luers/sample_model%v_p%v.txt", model, saveProp)).Done()
			if werr != nil {
			   fmt.Printf("Failed to write sample to disk\n%v", werr)
			}


		} else {
			fmt.Printf("Could not save a sample, number of observations not known")
		}
	}

	return da
}

type rec struct {
	concordance float64
	l2w         float64
	result      *duration.PHResults
	kr          *statmodel.KnockoffResult
}


func dropDrugs(da dstream.Dstream, dropProp float64) dstream.Dstream {
     fmt.Printf("Dropping drug therapeutic groups with less than %v proportion of positive values\n", dropProp)
     fmt.Printf("\nvariables: %v\n", da.Names())
     var droptg []string
     var dropvars []string
     tgnpos := make(map[string]int)
     da.Reset()
     ntot := 0
     fpos := func(r float64) bool {
     	  return r > 0.0
     }
     for da.Next() {
		for k, v := range da.Names() {
		       if strings.HasPrefix(v, "TG_") {
		       	  z := da.GetPos(k).([]float64)
			  ntot += len(z)
			  tgnpos[v] += floats.Count(fpos, z)
		       }
		}
     }
     
     
     for k, v := range tgnpos {
     	 tgprop := float64(v) / float64(ntot)
	 fmt.Printf("Drug: %v, Proportion of positive values: %v\n", k, tgprop)
     	 if tgprop < dropProp {
	    droptg = append(droptg, k)
	 }
     }
     fmt.Printf("Dropping main effects and interactions for these drugs with very few 1's:\n%v\n", droptg)
     
     if len(droptg) > 0 {
     	for _, v := range da.Names() {
	    for _, tgname := range droptg {
	    	// assuming no drug-drug interactions
	    	if strings.Contains(v, tgname) {
		   dropvars = append(dropvars, v)
		}
	    }
	    
	}
	fmt.Printf("dropvars = %v\n", dropvars)
	da = dstream.DropCols(da, dropvars...)
     }
     return da
}

func ridge(mode int, l2w []float64, ko bool, fl_save bool, fl_fullrank bool, fl_qr bool, dropDrugProp float64) error {

	da := modeldata(mode, fl_save, fl_fullrank, fl_qr, ko)
	
	if mode == 3 || mode == 4 || mode == 7 || mode == 8 {
	   da = dropDrugs(da, dropDrugProp)
	   da.Reset()
	   fmt.Printf("\nVariable names after dropping drug groups with low proportion of positive values: %v\n", da.Names())
	   fmt.Printf("---Checking weight and Age columns after dropping some drugs---\n")
	   chunknumber := 0
	   for da.Next() {
	       fmt.Printf("chunk number %v\n", chunknumber)
	       chunknumber++
	       aa := da.Get("Age").([]float64)
	       ww := da.Get("Weight").([]float64)
	       fmt.Printf("len(age) = %v, len(weight) = %v\n", len(aa), len(ww))
	   }
	   
	}

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
		//ko, err := statmodel.NewKnockoff(da, names)
		fmt.Printf("\nCreating Knockoff model with %v variables", len(names))
		fmt.Printf("\nNames of variables passed to Knockoff: %v\n", names)

		ko := statmodel.NewKnockoff(da, names)
		// if  err != nil {
		//	return err
		// }
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
			var dx dstream.Dstream
			// dx := dstream.Shallow(da)
			dx = da

			model := duration.NewPHReg(dx, "Time", "HF").OptSettings(opt).L2Weight(l2wgt).Weight("Weight")
			if !ko {
				// With knockoff we are already using normed data
				model = model.Norm()
			}
			fmt.Printf("Names inside PHReg: %v\n", model.DataSet().Names())
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
	tn := time.Now()
	ts := tn.Format("Jan2-15-04-05")
	fname_p := []string{fmt.Sprintf("coeff_%d_", mode), ts, ".txt"}
	fname := strings.Join(fname_p, "")
	if ko {
		fname = strings.Replace(fname, ".txt", "_ko.txt", 1)
	}
	if fl_fullrank {
		fname = strings.Replace(fname, ".txt", "_fullrank.txt", 1)
	}
	if fl_qr {
		fname = strings.Replace(fname, ".txt", "_qr.txt", 1)
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

	var ko, fl_fullrank, fl_qr, fl_save bool
	flag.BoolVar(&ko, "knockoff", false, "Use knockoff method")
	flag.BoolVar(&fl_save, "save", false, "Save sample of records")
	flag.BoolVar(&fl_fullrank, "fullrank", false, "Find maximal set of linearly independent columns")
	flag.BoolVar(&fl_qr, "qr", false, "Use rank-revealing QR to drop redundant columns")
	flag.Parse()
	fmt.Printf("ko: %v\nfullrank: %v\nqr: %v\nsave records: %v\n", ko, fl_fullrank, fl_qr, fl_save)
	data = dstream.NewBCols("data", 100000).Done()
	genvars()
	data = center(data)

	stats := dstream.Describe(data)
	for k, v := range stats {
		fmt.Printf("%-20s %f\n", k, v.SD)
	}
	dropDrugProp := 0.0001
	// l2w := []float64{0.1} 
	l2w := []float64{0.05, 0.1, 0.2, 0.4, 0.8, 1.6}
	//fmt.Printf("Restricting to a single hyperparameter value\n")
	for k := 4; k < 6; k++ {
		err := ridge(k, l2w, ko, fl_save, fl_fullrank, fl_qr, dropDrugProp)
		if err != nil {
			print(err)
		}
	}
}
