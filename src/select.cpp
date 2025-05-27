#include "select.hpp"

Select::Select(){}

Select::Select(Vect v, CritClust b, SelectReg sReg, int packSize)
{
    this->v = v;
    this->b = b;
    this->sReg = sReg;
    this->packSize = packSize;
}

Select::Select(Vect v, SelectReg sReg, int packSize)
{
    this->v = v;
    this->sReg = sReg;
    this->packSize = packSize;
}

List Select::selectS(std::vector<int> Order) {
    int InitialProjectsNb = this->v.experiments.size();
    double CritValue;
    const int numeromodeleaux = 1;
    int firstIndex = 1;
    int lastIndex = firstIndex + packSize;
    int ClustVar = 1000;

    std::vector<int> aux, varSelectReg, varSelectClust_aux, varSelectClust;
    double critClustaux = 0.0, critDiffClust = 0.0;
    List Mylist, Mylistaux;
    DataFrame missingVals;

    // Initialize variable selection with first order variable
    varSelectClust.push_back(Order[0]);
    Mylist = b.ClustBestModel(varSelectClust);

    // Initial criterion value and error check
    CritValue = as<double>(Mylist["criterionValue"]);
    std::string s_target = as<std::string>(Mylist["error"]);

    // Capture missing values from the initial model
    if (Mylist.containsElementNamed("missingValues")) {
        missingVals = as<DataFrame>(Mylist["missingValues"]);
    } else {
        // Create an empty DataFrame if no missing values
        missingVals = DataFrame::create(
            Named("row") = IntegerVector(0),
            Named("col") = IntegerVector(0),
            Named("value") = NumericVector(0)
        );
    }

    // Variable selection process
    while ((ClustVar > 0) && (firstIndex < (int)Order.size()) && (s_target == "No error")) {
        ClustVar = 0;
        for (int idx = firstIndex; idx < lastIndex; ++idx) {
            if (idx < (int)Order.size()) {
                // Reset temporary vectors
                aux.clear();
                varSelectReg.clear();
                varSelectClust_aux.clear();

                // Add current variable
                aux.push_back(Order[idx]);

                // Select regression variables
                varSelectReg = this->sReg.selectReg(varSelectClust, aux, InitialProjectsNb);

                // Augment variable selection
                varSelectClust_aux = this->v.ajouter_var(varSelectClust, aux);

                // Cluster with augmented variables
                Mylistaux = b.ClustBestModel(varSelectClust_aux);

                // Check model error
                s_target = as<std::string>(Mylistaux["error"]);
                if (s_target != "No error") continue;

                // Calculate BIC for regression
                List mylist = this->v.bicReggen(aux, varSelectReg, numeromodeleaux);

                // Compute criterion difference
                critClustaux = as<double>(Mylistaux["criterionValue"]);
                critDiffClust = critClustaux - CritValue - as<double>(mylist["bicvalue"]);

                // Update if improvement found
                if (critDiffClust > 0) {
                    varSelectClust = varSelectClust_aux;
                    Mylist = Mylistaux;
                    CritValue = critClustaux;

                    // Update missing values if available
                    if (Mylistaux.containsElementNamed("missingValues")) {
                        missingVals = as<DataFrame>(Mylistaux["missingValues"]);
                    } else {
                        // Create an empty DataFrame if no missing values
                        missingVals = DataFrame::create(
                            Named("row") = IntegerVector(0),
                            Named("col") = IntegerVector(0),
                            Named("value") = NumericVector(0)
                        );
                    }

                    ClustVar++;
                }
            }
        }

        // Move to next pack of variables
        firstIndex = lastIndex;
        lastIndex += packSize;
    }

    // Return results including missing values
    return List::create(
        Named("S") = wrap(varSelectClust),
        Named("model") = Mylist["model"],
        Named("criterionValue") = Mylist["criterionValue"],
        Named("criterion") = Mylist["criterion"],
        Named("nbcluster") = Mylist["nbcluster"],
        Named("parameters") = Mylist["parameters"],
        Named("proba") = Mylist["proba"],
        Named("partition") = Mylist["partition"],
        Named("missingValues") = missingVals
    );
}


vector<int> Select::selectW(vector<int> Order, vector<int>OtherVar)
{
    vector<int> varIndep;
    int InitialProjectsNb = this->v.experiments.size();
    int firstIndex, lastIndex, idx;
    lastIndex = Order.size();
    firstIndex = lastIndex - packSize;
    int ClustVar = 1000;
    
    //vector<int> OtherVar, varIndep_Role(Order.size()), varRegAux, aux;
    vector<int> varIndep_Role(Order.size()), varRegAux, aux;
    for(int l = 0; l < (int)Order.size(); ++l)
        varIndep_Role[l] = 0;
    
    varIndep.clear();
    while((ClustVar > 0) && (firstIndex >= 0))
    {
        ClustVar =0;
        for(idx = (lastIndex - 1); idx >= firstIndex; idx--)
        {
            if(idx >= 0)
            {
                varRegAux.clear(); aux.clear();
                //OtherVar.clear();
                //for(int l = 0; l < idx; ++l)
                //    if((varIndep_Role[l] == 0) && (l != idx))
                //        OtherVar.push_back(Order[l]);
                
                aux.push_back(Order[idx]);
                varRegAux=sReg.selectReg(OtherVar,aux,InitialProjectsNb);              	    
                if(varRegAux.empty())
                {
                    varIndep.push_back(Order[idx]);
                    ClustVar++;
                    varIndep_Role[idx] = 1; 
                } 
                
            } 
        } 
        lastIndex = firstIndex; 
        firstIndex -= packSize; 
    }
    return(varIndep);
};




