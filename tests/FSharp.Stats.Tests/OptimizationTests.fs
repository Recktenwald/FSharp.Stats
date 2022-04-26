module OptimizationTests

open Expecto

open FSharp.Stats
open FSharp.Stats.Optimization


[<Tests>]
let NelderMead = 
    testList "Nelder Mead" [
        testCase "Squared Distance from Origin" (fun () -> 
            let offset = 0.05
            let startingPoint = vector [1.0;2.0]
            let f v = Vector.dot v v
            let points = [|vector [1.0;2.0]; vector [1.0 + offset;2.0]; vector [1.0;2.0+2.0*offset]|]
            let values = Array.map f points
            let input = 
                Array.zip points values 
                |> Array.sortBy snd 
        
            let expectedMinimum = [|0.0;0.0|]
            let nmMinimum = 
                NelderMead.temp_doNSteps 200 1.0 2.0 0.5 0.5 f input
                |> fst
            Expect.floatClose Accuracy.low nmMinimum[0] expectedMinimum[0] "Minimums should be equal"
            Expect.floatClose Accuracy.low nmMinimum[0] expectedMinimum[1] "Minimums should be equal"
        )
    ]