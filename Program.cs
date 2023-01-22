// See https://aka.ms/new-console-template for more information
using ILOG.CPLEX;
using ILOG.Concert;

Console.WriteLine("Solve 2BPPDD");
getcplex();

static void getcplex()
{
    //Parameters
    double[] sizeb={10,10}; //W-H
    double[,] sizei={{9	, 5},
                     {2	, 4},
                     {6	,10},
                     {7	, 5},
                     {3	, 6},
                     {7	,10},
                     {5	, 1},
                     {5	, 3},
                     {9	, 6},
                     {4	, 2},
                     {7	, 6},
                     {2	, 7},
                     {3	, 8},
                     {10, 4},
                     {5	, 4},
                     {3	,10},
                     {3	, 8},
                     {8	, 7},
                     {3	, 8},
                     {7	, 8}};
    
    int n=sizei.GetLength(0);
    int nstar=2*n;
    double[,] sizeiprim=new double[nstar,sizei.GetLength(1)];
    for (int i = 0; i < n; i++)
    {
      sizeiprim[i,0]=sizei[i,0];
      sizeiprim[i,1]=sizei[i,1];
      sizeiprim[i+n,0]=sizei[i,1];
      sizeiprim[i+n,1]=sizei[i,0];
    }

    int bn=20; //number of the bins
    double [] di={418,395,276,279,419,305,291,238,283,149,
                  402,330,264,138,253,365,181,259,345,197};
    double [] diprim=new double[nstar];
    for (int i = 0; i < n; i++)
    {
        diprim[i]=di[i];
        diprim[i+n]=di[i];
    }
    double pb=100; //procceissing time of the bin


    Cplex mycplex=new Cplex();
    
    //Variables

    INumVar [] x=new INumVar[nstar];
    //IIntVar [] x=new IIntVar[nstar];
    INumVar [] y=new INumVar[nstar];
    //IIntVar [] y=new IIntVar[nstar];
    for (int i = 0; i < nstar; i++)
    {
        //x[i]=mycplex.NumVar(0,System.Double.MaxValue,"x"+"["+i+"]");
        x[i]=mycplex.NumVar(0,System.Double.MaxValue);
        //x[i]=mycplex.IntVar(0,10);
        //y[i]=mycplex.NumVar(0,System.Double.MaxValue,"y"+"["+i+"]");
        y[i]=mycplex.NumVar(0,System.Double.MaxValue);
        //y[i]=mycplex.IntVar(0,10);

    }

    INumVar [,] f=new INumVar[nstar,bn];
    for (int i = 0; i < nstar; i++)
    {
        for (int j = 0; j < bn; j++)
        {
            //f[i,j]=mycplex.BoolVar("f"+"["+i+","+j+"]");
            f[i,j]=mycplex.BoolVar();
        }
    }

    INumVar [,] l=new INumVar[nstar,nstar];
    INumVar [,] b=new INumVar[nstar,nstar];
    for (int i = 0; i < nstar; i++)
    {
        for (int j = 0; j < nstar; j++)
        {
            if ((i!=j) && (i!=(j+n)) && (j!=(i+n)))
            {
                //l[i,j]=mycplex.BoolVar("l"+"["+i+","+j+"]");
                l[i,j]=mycplex.BoolVar();
                //b[i,j]=mycplex.BoolVar("b"+"["+i+","+j+"]");
                b[i,j]=mycplex.BoolVar();
            }
        }
    }

    //INumVar lmax=mycplex.NumVar(System.Double.MinValue,System.Double.MaxValue,"lmax");
    INumVar lmax=mycplex.NumVar(System.Double.MinValue,System.Double.MaxValue);

    //Objective Function

    ILinearNumExpr obj=mycplex.LinearNumExpr();
    obj.AddTerm(lmax,1.0);
    mycplex.AddMinimize(obj);

    //Constraints
    //1
    for (int i = 0; i < nstar; i++)
    {
        for (int j = 0; j < nstar; j++)
        {
            if ((i<j) && (j!=(i+n)))
            {
                for (int k = 0; k < bn; k++)
                {
                    ILinearNumExpr expr=mycplex.LinearNumExpr();
                    expr.AddTerm(1.0,l[i,j]);
                    expr.AddTerm(1.0,l[j,i]);
                    expr.AddTerm(1.0,b[i,j]);
                    expr.AddTerm(1.0,b[j,i]);
                    expr.AddTerm(-1.0,f[i,k]);
                    expr.AddTerm(-1.0,f[j,k]);
                    mycplex.AddGe(expr,-1.0);
                }
                
            }
        }
    }

    //2
    for (int i = 0; i < nstar; i++)
    {
        for (int j = 0; j < nstar; j++)
        {
            if ((i!=j) && (i!=(j+n)) && (j!=(i+n)))
            {
                ILinearNumExpr expr=mycplex.LinearNumExpr();
                expr.AddTerm(1.0,x[i]);
                expr.AddTerm(-1.0,x[j]);
                expr.AddTerm(sizeb[0],l[i,j]);
                mycplex.AddLe(expr,(sizeb[0])-sizeiprim[i,0]);
            }
        }
    }

    //3
    for (int i = 0; i < nstar; i++)
    {
        for (int j = 0; j < nstar; j++)
        {
            if ((i!=j) && (i!=(j+n)) && (j!=(i+n)))
            {
                ILinearNumExpr expr=mycplex.LinearNumExpr();
                expr.AddTerm(1.0,y[i]);
                expr.AddTerm(-1.0,y[j]);
                expr.AddTerm(sizeb[1],b[i,j]);
                mycplex.AddLe(expr,(sizeb[1])-sizeiprim[i,1]);
            }
        }
    }

    //4
    for (int i = 0; i < nstar; i++)
    {
        ILinearNumExpr expr=mycplex.LinearNumExpr();
        expr.AddTerm(1.0,x[i]);
        mycplex.AddLe(expr,sizeb[0]-sizeiprim[i,0]);
    }

    //5
    for (int i = 0; i < nstar; i++)
    {
        ILinearNumExpr expr=mycplex.LinearNumExpr();
        expr.AddTerm(1.0,y[i]);
        mycplex.AddLe(expr,sizeb[1]-sizeiprim[i,1]);
    }

    //6
    for (int i = 0; i < n; i++)
    {
        ILinearNumExpr expr=mycplex.LinearNumExpr();
        for (int k = 0; k < bn; k++)
        {
            expr.AddTerm(1.0,f[i,k]);
            expr.AddTerm(1.0,f[i+n,k]);
        }
        mycplex.AddEq(expr,1.0);
    
    }

    //7
    for (int i = 0; i < nstar; i++)
    {
        ILinearNumExpr expr=mycplex.LinearNumExpr();
        for (int k = 0; k < bn; k++)
        {
            expr.AddTerm(((k+1)*pb)-diprim[i],f[i,k]);
        }
        expr.AddTerm(-1.0,lmax);
        mycplex.AddLe(expr,0);
    }
    

    
    //Solve
    if(mycplex.Solve())
    {
        System.Console.WriteLine("Objective Function={0}",mycplex.GetObjValue());

    }


}
