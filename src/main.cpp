#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <math.h>
#include <ctime>
#include <list>

#include "InstanceUCP.h"
#include "ModeleFlot.h"
#include "ModeleUCP.h"
#include "Process.h"

using namespace std ;




ILOSTLBEGIN

int process(InstanceProcessed I, ofstream & fichier, double & time, Methode met, IloEnv env) {

    // cout << "ici : " << met.getNum() << endl ;

    string nom = I.fileName() ;
    const char* file = nom.c_str() ;


    int id=I.id ;

    //IloEnv env ;
    InstanceUCP* inst = new InstanceUCP(env, file) ;
    cout << "nbG:" << inst->nbG << endl ;

    int n = inst->n ;
    int T = inst->T ;

    IloBoolVarArray x(env,n*T);
    IloBoolVarArray u(env,n*T);

    IloBoolVarArray f(env,n*(T+2)*(T+2));


    IntervalModel modelInt(env, inst);
    int Lmin = modelInt.Lmin;
    IloIntVarArray Y(env,inst->nbG*T*(T-Lmin+1), 0, n) ;

    IloModel model ;

    int ramp = met.Ramping();



    if (met.IneqVarY()) {
        if (!met.IneqSum()) {
            model = defineModel_y(env,inst,x,u) ;
        }
        else {
            model = defineModel_sum(env,inst, x,u, -5) ;
        }
    }
    else if (met.IneqSum()) {
        model = defineModel_sum(env,inst, x,u, -3) ;

    }

    else if (met.NumberOfOnes()) {
        model = defineModel_numberOfOnes(env,inst, x,u) ;
    }

    else if (met.AggregatedModel()) {
        model = AggregatedModel(env, inst) ;
    }

    else if (met.ModeleFlot()) {
        ModeleFlot flot(env, inst) ;
        model = flot.AggregatedFlowModel();
    }

    else if (met.ModeleIntervalle()) {
        model = modelInt.defineIntervalModel(Y) ;
    }

    else {
        model = defineModel(env,inst,x,u, 0, ramp) ;
        if (met.RSUonly()) {
            AddRSUIneq(model, env, inst,x,u,0);
        }
    }

    IloCplex cplex(model) ;


    //Paramètres
    cplex.setParam(IloCplex::Param::ClockType, 1); //1 : CPU TIME
    cplex.setParam(IloCplex::EpGap, 0.0000001) ;
    cplex.setParam(IloCplex::Param::TimeLimit, 3600) ;


    //Résolution et affichage de la solution

    cplex.solve();

    if (cplex.isPrimalFeasible()) {
        //cout << "feasible : " << cplex.isPrimalFeasible() << endl ;
        /*cplex.getValues(sub.x_frac, x) ;

    for (int t=0 ; t < T ; t++) {
        cout << "t=" << t << "    " ;
        for (int i= 0 ; i <sub.n ; i++) {
            cout << abs(sub.x_frac[i*T + t]) << " " ;
            if (sub.Last[i]) {
                cout << "         " ;
            }
        }
        cout << endl ;
    }
    cout << endl;*/

        ///Affichage solution optimale


        //    IloNumArray x_frac(env,0);
        //    cplex.getValues(x_frac, x) ;
        //    for (int i=0; i<n; i++) {
        //        cout << "unité " << i << " : " ;
        //        for (int t=0 ; t < T ; t++) {
        //            cout << x_frac[i*T +t] << " " ;
        //        }
        //        cout << endl ;
        //    }
        //    cout << endl ;

        //    IloNumArray u_frac(env,0);
        //    cplex.getValues(u_frac, u) ;
        //    for (int i=0; i<n; i++) {
        //        cout << "start up unité " << i << " : " ;
        //        for (int t=0 ; t < T ; t++) {
        //            cout << u_frac[i*T +t] << " " ;
        //        }
        //        cout << endl ;
        //    }
        //    cout << endl ;

        /* else {

    }*/

        double t = cplex.getCplexTime();
        double opt = cplex.getObjValue() ;

        fichier << met.getNum() <<  " & " << n << " & " << T  << " & " << I.symetrie << " & " << inst->nbG  << " & " << inst->MaxSize << " & " << inst->MeanSize  << " & " << id ;
        fichier << " & " << cplex.getObjValue()  ; //Optimal value
        fichier << " & " << cplex.getMIPRelativeGap() << " \\% " ; //approx gap
        fichier << " & " << cplex.getNnodes() ;
        fichier << " & " << t - time ;
        fichier <<" \\\\ " << endl ;

        time = t ;
    }


    //Destructeurs
    // delete inst ;
    // delete dataNode ;
    // env.end() ;

    return 1;
}


int
main(int argc,char**argv)
{
    srand(time(NULL));

    //définition des méthodes de résolution
    Methode DefaultCplex ;

    Methode IneqPures;
    IneqPures.UseIneqSum();

    Methode RampModel;
    RampModel.UseRampConstraints();

    Methode RampIneqRSU;
    RampIneqRSU.UseRampConstraints();
    RampIneqRSU.AddIneqRSU() ;

    Methode ModeleIntervalle ;
    ModeleIntervalle.UseModeleInterval();

    Methode IneqVarY;
    IneqVarY.UseIneqVarY();

    Methode Flot;
    Flot.UseModeleFlot();

    Methode IneqNumberOfOnes;
    IneqNumberOfOnes.UseNumberOfOnes();

    Methode IneqSumAndVarY;
    IneqSumAndVarY.UseIneqVarY();
    IneqSumAndVarY.UseIneqSum();
    IneqSumAndVarY.setNum(-5);

    Methode AggregModel;
    AggregModel.UseAggregatedModel();

    /////////////////////// SI ARGUMENTS //////////////////////
    if (argc>1) {
        ofstream fichier("result.txt", std::ofstream::out | std::ofstream::app);

        int met= atoi(argv[1]);
        string localisation = argv[2] ;

        int n = atoi(argv[3]);
        int T = atoi(argv[4]);
        int bloc = atoi(argv[5]);
        int demande = atoi(argv[6]);
        int sym = atoi(argv[7]);

        int cat01 = atoi(argv[8]);
        int intra = atoi(argv[9]);

        int id = atoi(argv[10]);

        InstanceProcessed Instance = InstanceProcessed(n, T, bloc, demande, sym, cat01, intra, id, localisation) ;

        double time = 0 ;
        IloEnv env ;

        if (met==1) {
            env=IloEnv() ;
            process(Instance, fichier, time, RampModel , env) ;
            env.end() ;

            env=IloEnv() ;
            process(Instance, fichier, time, RampIneqRSU, env) ;
            env.end() ;

            env=IloEnv() ;
            process(Instance, fichier, time, ModeleIntervalle, env) ;
            env.end() ;
            fichier << endl ;
        }
    }


    ////////////////////////// SI PAS D'ARGUMENTS ////////////////////
    if (argc==1) {
        ofstream fichier("result.txt");

        //fichier << "Instance & n & T & Sym & nG & max & mean & OptVal & RootRelax & ApproxGap &  USCuts & IUSCuts & SepTime & Prof & Nodes & CplexCuts & CPU \\\\ " << endl;
        fichier << "Methode & n & T & OptVal & Nodes & NbFixs & TimeFix & CPU \\\\ " << endl;

        double time = 0 ;

        //Paramètres de l'instance

        int T = 96;
        int n = 60 ;
        int sym = 2 ;
        int demande = 3;
        int cat01 = 0;
        int bloc = 1;

        int intra =0 ;

        string localisation = "data/smaller/" ;
        InstanceProcessed Instance = InstanceProcessed(n, T, bloc, demande, sym, cat01, intra, 0, localisation) ;

        fichier << localisation << endl ;
        Instance.localisation = localisation ;

        n=20;
        T=24;
        Instance.n=n;
        Instance.T=T ;
        IloEnv env ;

        for (sym= 3; sym >=3 ; sym--) {
            Instance.symetrie = sym ;

            for (T=24 ; T <= 965 ; T*=2) {

                Instance.T=T ;


                for (int id=1; id <=20; id++) {
                    Instance.id = id ;

                    env=IloEnv() ;
                    cout <<"start ramp model" << endl ;
                    process(Instance, fichier, time, RampModel , env) ;
                    cout <<"end ramp model" << endl ;
                    env.end() ;

                    env=IloEnv() ;
                    process(Instance, fichier, time, RampIneqRSU, env) ;
                    env.end() ;

                    env=IloEnv() ;
                    process(Instance, fichier, time, ModeleIntervalle, env) ;
                    env.end() ;

                    fichier << endl ;
                }

                fichier << endl ;
                fichier << endl ;
            }
        }
    }


    //        int intra = 0 ;
    //        string localisation = "data/Litt_Real/" ;
    //        InstanceProcessed Instance = InstanceProcessed(n, T, bloc, demande, sym, cat01, intra, 0, localisation) ;

    //        fichier << localisation << endl ;
    //        Instance.localisation = localisation ;

    //        n=30;
    //        T=96;
    //        Instance.n=n;
    //        Instance.T=T ;
    //        IloEnv env ;

    //        for (sym= 3; sym >=3 ; sym--) {
    //            Instance.symetrie = sym ;

    //            for (int id=1; id <=10; id++) {
    //                Instance.id = id ;

    //                env=IloEnv() ;
    //                cout <<"start ramp model" << endl ;
    //                process(Instance, fichier, time, RampModel , env) ;
    //                cout <<"end ramp model" << endl ;
    //                env.end() ;

    //                env=IloEnv() ;
    //                process(Instance, fichier, time, RampIneqRSU, env) ;
    //                env.end() ;

    //                /*env=IloEnv() ;
    //            process(Instance, fichier, time, ModeleIntervalle, env) ;
    //            env.end() ;*/

    //                fichier << endl ;
    //            }

    //            fichier << endl ;
    //            fichier << endl ;
    //        }
    //    }




    return 0 ;
}
