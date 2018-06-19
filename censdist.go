package main

import (
	"fmt"
	"os"

	"github.com/kshedden/dstream/dstream"
	"github.com/kshedden/duration"
)

const (
	minage float64 = 50
)

func main() {

	data := dstream.NewBCols("data", 1000000).Done()

	agefilter := func(x interface{}, keep []bool) bool {
		age := x.([]float64)
		for i, a := range age {
			if a < minage {
				keep[i] = false
			}
		}
		return true
	}
	data = dstream.Filter(data, map[string]dstream.FilterFunc{"Age": agefilter})

	rev := func(v map[string]interface{}, x interface{}) {
		status := v["HF"].([]float64)
		z := x.([]float64)

		for i := range z {
			z[i] = 1 - status[i]
		}
	}

	data = dstream.Apply(data, "Rstatus", rev, "float64")
	data.Reset()
	data = dstream.MemCopy(data)

	sf := duration.NewSurvfuncRight(data, "Time", "Rstatus").Done()

	ti := sf.Time()
	sp := sf.SurvProb()

	fid, err := os.Create(fmt.Sprintf("censdist_%.0f.csv", minage))
	if err != nil {
		panic(err)
	}
	defer fid.Close()

	for i := range ti {
		fid.Write([]byte(fmt.Sprintf("%f,%f\n", ti[i], sp[i])))
	}
}
