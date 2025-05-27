#include "critClust.hpp"

//****************************************************************************//
// Constructor ***************************************************************//
CritClust::CritClust() {}

CritClust::CritClust(int k, std::string framework, std::string model_name, NumericMatrix data, 
                     std::string crit, IntegerVector knownlabels, bool DA)
{
    this->crit = crit;
    this->framework = framework;
    this->model_name = model_name;
    this->k = k;
    this->data = data;
    this->knownlabels = knownlabels;
    this->DA = DA;
}

List CritClust::ClustBestModel(std::vector<int> numExp)
{
    // Default error result to return if needed
    List defaultResult = List::create(
        Named("criterionValue") = NA_REAL,
        Named("criterion") = crit,
        Named("nbcluster") = k,
        Named("model") = model_name,
        Named("parameters") = R_NilValue,
        Named("proba") = R_NilValue,
        Named("partition") = R_NilValue,
        Named("error") = "Model fitting failed",
        Named("missingValues") = DataFrame::create(
            Named("row") = NumericVector(),
            Named("col") = NumericVector(),
            Named("value") = NumericVector()
        ),
        Named("S") = NumericVector() 
    );

    try {
        Environment base("package:base");
        Function dataframe = base["data.frame"];
        Function parseR = base["parse"];
        Function evalR = base["eval"];
        
        // Prepare subset of data
        NumericMatrix dataAux(data.nrow(), numExp.size());
        for (size_t j = 0; j < numExp.size(); ++j) {
            dataAux(_, j) = data(_, numExp[j] - 1);
        }
        DataFrame df = dataframe(dataAux);

        // Framework decision 
        if (framework == "Mclust") {
            if (DA) {
                Rcout << "Warning: Mclust does not support supervised learning. Using unsupervised.\n";
            }
            Environment mclustEnv("package:mclust");
            Function MclustFunc = mclustEnv["Mclust"];
            Function hcFunction = mclustEnv["hc"];

            Rcpp::DataFrame Rdf = Rcpp::as<Rcpp::DataFrame>(df); 
            int numCols = Rdf.ncol();

            // CharacterVector hcModelForAgglomeration;
            // List initializationList; 
            

            // if (numCols == 1) {                    
            //     bool hasV = (model_name.find('V') != std::string::npos);
            //     hcModelForAgglomeration = CharacterVector::create(hasV ? "V" : "E");
            // } else {                               
            //     hcModelForAgglomeration = CharacterVector::create("VVV");
            // }
                 
            // try {
            //     List hcResult = hcFunction(
            //         Named("data")      = df,
            //         Named("modelName") = hcModelForAgglomeration);
            //     initializationList = List::create(Named("hcPairs") = hcResult);
            // } catch (...) {
            //     initializationList = List::create();
            // }

            List mclustArgs;
            mclustArgs["data"] = df;
            mclustArgs["G"] = k;
            if (numCols == 1) {
                 mclustArgs["modelNames"] = CharacterVector::create("E", "V");
            } else {
                mclustArgs["modelNames"] = model_name;
            }
            // if (initializationList.size() > 0  && initializationList.containsElementNamed("hcPairs"))
            //     mclustArgs["initialization"] = initializationList;
            mclustArgs["verbose"] = false;
            
            List mclustResult;
            try {
                Function do_call = base["do.call"];
                mclustResult = do_call(MclustFunc, mclustArgs);

            } catch (std::exception &ex) {
                std::string error_msg = std::string("Mclust error: ") + ex.what();
                Rcout << error_msg << "\n";
                defaultResult["error"] = error_msg;
                return defaultResult;
            } catch (...) {
                std::string error_msg = "Unknown error during Mclust call.";
                Rcout << error_msg << "\n";
                defaultResult["error"] = error_msg;
                return defaultResult;
            }

            DataFrame missingVals = DataFrame::create(
                Named("row") = NumericVector(),
                Named("col") = NumericVector(), 
                Named("value") = NumericVector()
            );                                            
            return List::create(
                Named("criterionValue") = as<double>(mclustResult["bic"]),
                Named("criterion") = "BIC",
                Named("nbcluster") = as<int>(mclustResult["G"]),
                Named("model") = as<std::string>(mclustResult["modelName"]),
                Named("parameters") = mclustResult["parameters"],
                Named("proba") = mclustResult["z"],
                Named("partition") = mclustResult["classification"],
                Named("error") = "No error",
                Named("missingValues") = missingVals
            );
        }
        // ------------------ Rmixmod Branch --------------------- //
        else if(framework == "Rmixmod") {
            Environment Rmixmod("package:Rmixmod");
            Function RmixmodLearn = Rmixmod["mixmodLearn"];
            Function RmixmodCluster = Rmixmod["mixmodCluster"];
            Function RmixmodStrategy = Rmixmod["mixmodStrategy"];
            Function RmixmodGaussianModel = Rmixmod["mixmodGaussianModel"];
            
            // Parse the model string to create a model object
            SEXP modelObject;
            try {
                SEXP parsedModel = parseR(Named("text") = model_name);
                modelObject = evalR(parsedModel);
            } catch(std::exception &ex) {
                modelObject = RmixmodGaussianModel(Named("family") = "all");
                Rcout << "Warning: Failed to parse model string '" << model_name << "': " << ex.what() << "\n";
                Rcout << "Using default all model instead.\n";
            } catch(...) {
                modelObject = RmixmodGaussianModel(Named("family") = "all");
                Rcout << "Warning: Failed to parse model string '" << model_name << "', using default all model.\n";
            }
            
            if(!DA) {
                // Unsupervised Clustering with Rmixmod
                S4 mixmodstrategy = RmixmodStrategy(
                    Named("nbTry") = 4,
                    Named("nbIterationInAlgo") = 200,
                    Named("nbTryInInit") = 100,
                    Named("nbIterationInInit") = 10,
                    Named("epsilonInAlgo") = 1e-8,
                    Named("initMethod") = "SEMMax"
                );
                
                CharacterVector criterionVec;
                if (crit == "BIC" || crit == "ICL"){
                    criterionVec = CharacterVector::create(crit);
                } else {
                    criterionVec = CharacterVector::create("BIC");
                    Rcout << "Warning: Unknown criterion '" << crit << "', using BIC instead.\n";
                }

                S4 xem = RmixmodCluster(
                    Named("data") = df,
                    Named("nbCluster") = k,
                    Named("models") = modelObject,
                    Named("strategy") = mixmodstrategy,
                    Named("criterion") = criterionVec
                );
                
                S4 bestResult = xem.slot("bestResult");
                DataFrame missingVals = DataFrame::create(
                    Named("row") = NumericVector(),
                    Named("col") = NumericVector(),
                    Named("value") = NumericVector()
                );
                
                std::string errorMsg = "No error";
                if (bestResult.hasSlot("error")) {
                    SEXP errorSlot = bestResult.slot("error");
                    if (!Rf_isNull(errorSlot)) {
                        errorMsg = as<std::string>(errorSlot);
                    }
                }
                
                List result = List::create(
                    Named("criterionValue") = -as<double>(bestResult.slot("criterionValue")),
                    Named("criterion") = bestResult.slot("criterion"),
                    Named("nbcluster") = bestResult.slot("nbCluster"),
                    Named("model") = bestResult.slot("model"),
                    Named("parameters") = bestResult.slot("parameters"),
                    Named("proba") = bestResult.slot("proba"),
                    Named("partition") = bestResult.slot("partition"),
                    Named("error") = errorMsg,
                    Named("missingValues") = missingVals
                );
                return result;
            }
            else {
                // Supervised learning with Rmixmod
                StringVector criteria;
                if (crit == "BIC" || crit == "CV") {
                    criteria = StringVector::create(crit);
                } else {
                    criteria = StringVector::create("BIC"); 
                    Rcout << "Warning: Using BIC criterion for supervised learning.\n";
                }
                
                S4 xem = RmixmodLearn(
                    Named("data") = df,
                    Named("knownLabels") = knownlabels,
                    Named("models") = modelObject,
                    Named("criterion") = criteria
                );
                
                S4 bestResult = xem.slot("bestResult");
                NumericVector critvalues = bestResult.slot("criterionValue");
                DataFrame missingVals = DataFrame::create(
                    Named("row") = NumericVector(),
                    Named("col") = NumericVector(),
                    Named("value") = NumericVector()
                );
                
                std::string errorMsg = "No error";
                if (bestResult.hasSlot("error")) {
                    SEXP errorSlot = bestResult.slot("error");
                    if (!Rf_isNull(errorSlot)) {
                        errorMsg = as<std::string>(errorSlot);
                    }
                }
                
                List result = List::create(
                    Named("criterionValue") = -critvalues[0],
                    Named("criterion") = bestResult.slot("criterion"),
                    Named("nbcluster") = bestResult.slot("nbCluster"),
                    Named("model") = bestResult.slot("model"), 
                    Named("parameters") = bestResult.slot("parameters"),
                    Named("proba") = R_NilValue, 
                    Named("partition") = bestResult.slot("partition"),
                    Named("error") = errorMsg,
                    Named("missingValues") = missingVals
                );
                return result;
            }
        }
        // ------------------ MixAll Branch (default) --------------------- //
        else {
            Environment MixAll("package:MixAll");
            Function clusterDiagGaussian = MixAll["clusterDiagGaussian"];
            Function clusterStrategy = MixAll["clusterStrategy"];
            Function learnDiagGaussian = MixAll["learnDiagGaussian"];
            Function missingValues_func = MixAll["missingValues"];
            
            S4 quick_precise_strategy = clusterStrategy(
                Named("nbTry") = 1,
                Named("nbInit") = 100,
                Named("initMethod") = "random",
                Named("initAlgo") = "SEM",
                Named("nbInitIteration") = 10,
                Named("initEpsilon") = 1e-4,
                Named("nbShortRun") = 10,
                Named("shortRunAlgo") = "EM",
                Named("nbShortIteration") = 100,
                Named("shortEpsilon") = 1e-7,
                Named("longRunAlgo") = "EM",
                Named("nbLongIteration") = 200,
                Named("longEpsilon") = 1e-8
            );
        
            S4 xem;
            CharacterVector models = CharacterVector::create(model_name);
            
            if(!DA){
                xem = clusterDiagGaussian(
                    Named("data") = df,
                    Named("nbCluster") = k,
                    Named("models") = model_name,
                    Named("strategy") = quick_precise_strategy, 
                    Named("criterion") = crit,
                    Named("nbCore") = 1
                );
            }
            else{
                IntegerVector labels = knownlabels;
                xem = learnDiagGaussian(
                    Named("data") = df,
                    Named("labels") = labels,
                    Named("models") = model_name,
                    Named("algo") = "simul",
                    Named("nbIter") = 100,
                    Named("epsilon") = 1e-08,
                    Named("criterion") = crit,
                    Named("nbCore") = 1
                );
            }
            
            if (xem.hasSlot("error")) {
                SEXP errorSlot = xem.slot("error");
                if (!Rf_isNull(errorSlot)) {
                    return List::create(Named("error") = as<std::string>(errorSlot));
                }
            }
            
            NumericMatrix imputedData = missingValues_func(xem);
            DataFrame missingVals = DataFrame::create(
                                        Named("row") = imputedData(_, 0),
                                        Named("col") = imputedData(_, 1),
                                        Named("value") = imputedData(_, 2)
                                    );
            
            List result = List::create(
                Named("criterionValue") = -as<double>(xem.slot("criterion")),
                Named("criterion") = xem.slot("criterionName"),
                Named("nbcluster") = xem.slot("nbCluster"),
                Named("model") = model_name,
                Named("parameters") = xem.slot("component"),
                Named("proba") = xem.slot("tik"),
                Named("partition") = xem.slot("zi"),
                Named("error") = "No error",
                Named("missingValues") = missingVals
            );
            return result;
        }
    }
    catch (std::exception &ex) {
        // return List::create(Named("error") = ex.what());
        defaultResult["error"] = ex.what();
        return defaultResult;
    }
    catch (...) {
        // return List::create(Named("error") = "Unknown error occurred in ClustBestModel");
        defaultResult["error"] = "Unknown error occurred in ClustBestModel";
        return defaultResult;
    }
}

