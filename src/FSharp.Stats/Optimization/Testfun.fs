namespace FSharp.Stats.Optimization

module Testfun = 
    open FSharp.Stats
    open System
    // Unbounded Optimization
    // // Rosenbrock function
    let rosenbrock v = 
        let n = Vector.length v 
        [0..(n-2)]
        |> List.map (fun i -> 100.0*( (v[i+1] - v[i]*v[i])**2.0) + ((1.0 - v[i])**2.0))
        |> List.sum 

    // // Sphere function
    let sphere v = 
        Vector.dot v v

    // Bounded Optimization 

    // // Ackley function 
    let ackley v = 
        let n = float (Vector.length v)
        let firstSummand = 
            20.0*exp(-0.2*sqrt((sphere v)/n))

        let sumOfCos = 
            v
            |> Vector.map (fun vi -> cos(2.0*Math.PI*vi))
            |> Vector.sum 
        let secondSummand = 
            exp(sumOfCos/n)
        
        - firstSummand - secondSummand + 20.0 + exp(1.0)

    // //  Himmelblau function
    /// This is a 2D -> 1D function.
    /// Only uses the first two entries if a longer vector is passed as an argument.
    let himmelblau (v:Vector<float>) = 
        (v[0]**2.0 + v[1] - 11.0)**2.0 + (v[0] + v[1]**2.0 - 7.0)**2.0
    