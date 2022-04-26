namespace FSharp.Stats.Optimization

module NelderMead =
    // Implementation of the Nelder-Mead Algorithm with parameters from 
    // Implementing the Nelder-Mead simplex algorithm with adaptive parameters by Gao and Han (2012)

    // In the n-dimensional problem, an n+1 dimensional simplex is used
    // and when sorting we fix the convention that f(x[0]) is the smallest value
    // and that the algorithm is used to minimize the function value.
    
    // We terminate the algorithm and return the best value so far when one of the following three conditions is satisfied
    // 1.
    //      a. All the function values are very close to the best one, i.e. max |f0-fi| < tolFun; and
    //      b. All the points are very close to the best one, i.e. max |x0-xi| < tolX (where we use the supremum metric)
    // 2. Number of iterations exceeds maxIter
    // 3. Number of function calls exceed maxFunEval


    open FSharp.Stats 
    
    type Vertex = Vector<float> * float

    /// Infinity distance of two vectors
    let inline supremum (v1:Vector<float>) (v2:Vector<float>) = 
        (v1-v2).InternalValues |> Array.map System.Math.Abs |> Array.max 

    let shouldTerminate tolFunOp tolXOp (simplex: Vertex array) : bool = 
        let xs, fs = Array.unzip simplex 
        let terminateFun = 
            match tolFunOp with 
            | None -> false
            | Some tolFun ->
                Array.tail fs 
                |> Array.map (fun fi -> fi - fs[0])
                |> Array.max
                |> fun x -> x < tolFun
        let terminatePoint = 
            match tolXOp with 
            | None -> false
            | Some tolX -> 
                Array.tail xs 
                |> Array.map (fun xi -> supremum xi xs[0])
                |> Array.max 
                |> fun x -> x < tolX
        terminateFun && terminatePoint

    /// Function that corresponds to one iteration
    /// The array that is given as an argument gets changed in place
    /// It is also assumed to be ordered already, 
    /// from the lowest function value at index 0 to the highest
    /// 
    /// For later use it also counts how many times the function has been called and returns that
    let private oneStepInPlace alpha beta gamma delta (f: Vector<float> -> float) (simplex_ : Vertex array) = 
        let mutable funEvals = 0
        
        let n = Array.length simplex_ - 1 
        let midPoint = 
            Vector.scale (1.0/(float (n))) 
                (Array.fold (fun state t -> state + (fst t)) (Vector.zeroCreate n) (Array.take n simplex_))
        let xr = midPoint + Vector.scale alpha (midPoint - (fst simplex_[n]))
        let fr = f xr 
        funEvals <- funEvals + 1
        let mutable doShrink = false 
        if fr < snd simplex_[0] then 
            // Expansion
            let xe = midPoint + Vector.scale beta (xr - midPoint)
            let fe = f xe
            funEvals <- funEvals + 1
            if fe < fr then 
                simplex_[n] <- (xe,fe)
            else 
                simplex_[n] <- (xr,fr)
        else if fr < snd simplex_[n-1] then 
            simplex_[n] <- (xr,fr)
        else if fr < snd simplex_[n] then 
            // Outside Contraction
            let xo = midPoint + Vector.scale gamma (xr - midPoint)
            let fo = f xo 
            funEvals <- funEvals + 1
            if fo < fr then 
                simplex_[n] <- (xo,fo)
            else   
                doShrink <- true
        else
            // Inside Contraction 
            let xi = midPoint - Vector.scale gamma (xr - midPoint) 
            let fi = f xi 
            funEvals <- funEvals + 1
            if fi < fr then 
                simplex_[n] <- (xi,fi)
            else
                doShrink <- true
        if doShrink then 
            let x0 = Array.head simplex_ |> fst
            // note that we shrink towards the best point so far, so x0 stays the same and we do n iterations, not n+1
            for i in [1..n] do 
                let xi = x0 + Vector.scale delta ((fst simplex_[i]) - x0) 
                simplex_[i] <- (xi, f xi)
            funEvals <- funEvals + n 
        else 
            ()
        // The obvious performance optimization is to insert the new element
        // in linear time, and only do a full sort after the shrinking step
        // which is (empirically) relatively rare.
        Array.sortInPlaceBy snd simplex_
        // return the number of function calls made
        funEvals

    let temp_doNSteps N alpha beta gamma delta f simplex_ =
        for i in [1..N] do 
            oneStepInPlace alpha beta gamma delta f simplex_ 
            |> ignore
        Array.head simplex_

    // Now that we have the base logic defined it is time to make a nicer API
    // In my opinion it improves readability to have the logic on top
    // and define the helper / config types afterwards. 
    
    /// tolFun : When |f0 - fn| < tolFun then terminate. Set to None if this should not be a criterion
    /// tolX : When max |x0-xi| < tolX then terminate. Set to None if this should not be a criterion
    /// maxIter : Do at most this many iterations.
    /// maxFunEval : Do at most this many function calls.
    ///     If maxIter and maxFunEval are both None, they will both be set to 200 * (dimension of the domain)
    ///     If exactly one of maxIter or maxFunEval is None, it will not be used to decide termination.
    /// alpha, beta, gamma, delta : Parameters in the Nelder-Mead algorithm. 
    ///     If set to None, they will be chosen depending on the dimension of the domain
    type Parameters = {
        tolFun : float option 
        tolX : float option 
        maxIter : int option 
        maxFunEval : int option 
        alpha : float option 
        beta : float option 
        gamma : float option 
        delta : float option 
    }

    let defaultParameters : Parameters = {
        tolFun = Some 1e-4
        tolX = Some 1e-4
        maxIter = None
        maxFunEval = None 
        alpha = None 
        beta = None 
        gamma = None 
        delta = None
    }

    
    // Adding some additional generic type parameters to mark this function as parametric
    let minimizeNelderMead<'a> (parameters : Parameters) (f : Vector<float> -> float) (initialGuess:'a) =
        // We need to count the function evaluations
        // Since we might need to call the function to construct the intial simplex, 
        // we define this at the very start
        let mutable funEvalCount = 0 
        
        // It is usually not recommended to match on types, but I think here it makes sense 
        // since we can allow for three different types of input
        // while keeping the api clean
        let dim, simplex = 
            match box initialGuess with 
            | :? array<Vertex> as initialSimplex -> 
                // In this case we are given a simplex.
                // We are checking that it has dimension + 1 points
                let dim = initialSimplex[0] |> fst |> Vector.length
                if dim + 1 <> Array.length initialSimplex then 
                    raise (System.ArgumentException("Initial simplex must consist of n+1 points, where n is the dimension of the domain."))
                dim, initialSimplex
            | :? array<Vector<float>> as initialPoints ->
                // In this case we are only given initial points and have to get the function values
                let dim = initialPoints[0] |> Vector.length
                if dim + 1 <> Array.length initialPoints then 
                    raise (System.ArgumentException("Initial simplex must consist of n+1 points, where n is the dimension of the domain."))
                funEvalCount <- funEvalCount + dim + 1
                let initialSimplex = 
                    initialPoints
                    |> Array.map (fun x -> (x, f x))
                dim, initialSimplex
            | :? Vector<float> as initialPoint ->
                // in this case we are only given one point and construct the simplex ourselves
                // If the i-th entry of the initial Point is 0 we add 0.00025 to that value, otherwise we add 0.05
                // To do this we need a helper function that creates the canonical basis vectors 
                // And since we are going to scale them anyway we just do this in the function 
                let dim = initialPoint |> Vector.length
                let basisVector dim i value = 
                    Vector.init dim (fun j -> if i=j then value else 0)
                
                // This fold looks at each entry in turn, remembering its index,
                // checks if it is close to zero or not, and adds a basis vector as described above
                // It then adds this to a running list of points
                let initialPoints = 
                    initialPoint 
                    |> Vector.foldi 
                        (fun index state xi ->
                            if System.Math.Abs xi < 1e-4 then 
                                (initialPoint + (basisVector dim index 0.00025))::state
                            else 
                                (initialPoint + (basisVector dim index 0.05))::state
                        ) [initialPoint]
                funEvalCount <- funEvalCount + dim + 1
                let initialSimplex = 
                    initialPoints
                    |> Array.ofList
                    |> Array.map (fun x -> (x, f x))
                dim, initialSimplex        
            | _ -> raise (System.ArgumentException("""The initial guess could not be interpreted. 
            Use an array of vectors, a single vector, or an array of tuples (vector, function value at that point)."""))

        

        // Set values for termination
        let maxIter, maxFunEval = 
            match parameters.maxIter, parameters.maxFunEval with 
            | None, None -> dim * 200, dim * 200
            | None, Some x -> System.Int32.MaxValue, x
            | Some x, None -> x, System.Int32.MaxValue
            | Some x, Some y -> x, y
        let alpha = defaultArg parameters.alpha 1.0 
        let beta = defaultArg parameters.beta (if dim > 1 then 2.0 else 1.0 + 2.0/(float dim))
        let gamma = defaultArg parameters.gamma (if dim > 1 then 0.5 else 0.75 - 1.0/(2.0*(float dim)))
        let delta = defaultArg parameters.delta (if dim > 1 then 0.5 else 1.0 - 1.0/(float dim))

        // sort and copy the simplex, so that we can do everything else in place, while keeping an immutable surface
        let simplex_ = Array.sortBy snd simplex 
        let mutable iterCount = 0 
        let mutable doTerminate = false 
        while iterCount < maxIter && not doTerminate do 
            funEvalCount <- funEvalCount + oneStepInPlace alpha beta gamma delta f simplex_
            doTerminate <- shouldTerminate parameters.tolFun parameters.tolX simplex_

        simplex_
        
    