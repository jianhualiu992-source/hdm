namespace hengduanmian;
public class MyFunctions : IExcelAddIn
{
    #region 填挖面积
    [ExcelFunction(Description = "填挖面积",Category ="横断面相关")]
    public static double[] hdm_CutAndFill(

       [ExcelArgument(Name = "sjx", Description = "设计点x，横向或纵向连续单元格")] double[] xsA,
       [ExcelArgument(Name = "sjy", Description = "设计点y，一维Range")] double[] ysA,
       [ExcelArgument(Name = "dmx", Description = "地面点x，一维Range")] double[] xsB,
       [ExcelArgument(Name = "dmy", Description = "地面点y，一维Range")] double[] ysB,
       [ExcelArgument(Name = "延长长度", Description = "比如1")] double epsilon)
    {
        return hdm.cutAndFillArea(xsA, ysA, xsB, ysB, epsilon,60);
    }
    #endregion

    #region k2xy
    [ExcelFunction(Description = "已知里程，计算对应的x y", Category = "横断面相关")]
    public static double[] hdm_k2xy(
        [ExcelArgument(Name = "name", Description = "\n第二段\n卡科拉")] string name,
        [ExcelArgument(Name = "k", Description = "桩号")] double k,
        [ExcelArgument(Name = "b", Description = "宽度，左负右正")] double b,
        [ExcelArgument(Name = "z", Description = "右夹角")] double z
        )
    {
        var pqx = Data.pqx[name];
        //System.Windows.Forms.MessageBox(pqx[0][1] + "");
        return hdm.Dantiaoxianludange(pqx, k, b, z);
    }
    #endregion

    #region xy2k
    [ExcelFunction(Description = "已知x,y计算对应的k b", Category = "横断面相关")]
    public static double[] hdm_xy2k(
        [ExcelArgument(Name = "name", Description = "\n第二段\n卡科拉")] string name,
        [ExcelArgument(Name = "x", Description = "x")] double x,
        [ExcelArgument(Name = "y", Description = "y")] double y)
    {
        var pqx = Data.pqx[name];
        //System.Windows.Forms.MessageBox(pqx[0][1] + "");
        return hdm.Fs(pqx, x,y);
    }
    #endregion

    #region h
    [ExcelFunction(Description = "已知里程，计算对应的h", Category = "横断面相关")]
    public static double hdm_center_h(  
        [ExcelArgument(Name = "k", Description = "桩号")] double k )
    {
        var sqxb = Data.sqxb;
        //System.Windows.Forms.MessageBox(pqx[0][1] + "");
        return Math.Round (hdm.h(k, sqxb),3);
    }
    #endregion


    #region slope
    [ExcelFunction(Description = "已知里程，计算对应的 横坡 ", Category = "横断面相关")]
    public static double[] hdm_slope(
        [ExcelArgument(Name = "k", Description = "桩号")] double k)
    {
        var slope = Data.slope;
        for (var i = 0; i < slope.Length - 1; i++)
        {
            if (k >= slope[i][0] && k <= slope[i + 1][0])
            {
                double L = slope[i][1] + (slope[i + 1][1] - slope[i][1]) / (slope[i + 1][0] - slope[i][0]) * (k - slope[i][0]);
                double R = slope[i][2] + (slope[i + 1][2] - slope[i][2]) / (slope[i + 1][0] - slope[i][0]) * (k - slope[i][0]);
                return [L, R];
            }
        }
        return [0.0, 0.0];
    }

    #endregion

    #region 结构层宽度
    [ExcelFunction(Description = "结构层宽度", Category = "横断面相关")]
    public static double[] hdm_结构层宽度(

       [ExcelArgument(Name = "k", Description = "里程")] double k,
       [ExcelArgument(Name = "kd", Description = "顶层宽度(左负右正)")] double kd,
       [ExcelArgument(Name = "hd", Description = "距离顶层厚度")] double hd)
    {
        var hp = hdm_slope(k);
        var hpL = hp[0];
        var hpR = hp[1];

        var h=hdm_center_h(k);
        if (kd < 0)
        {
            var y1 = h - kd * hpL;
            var y2 = h - hd;
            var jd1 = hdm.Intersections_line(kd, y1, kd - 2, y1 - 1,0,y2,-1,y2+hpL);
            return jd1;
        }
        else
        {
            var y1 = h + kd * hpR;
            var y2 = h - hd;
            var jd1 = hdm.Intersections_line(kd, y1, kd + 2, y1 - 1, 0, y2, 1, y2 + hpR);
            return jd1;
        }
    }
    #endregion

    #region 已知折线,点坐标x,计算对应的y
    [ExcelFunction(Description = "填挖判断", Category = "横断面相关")]
    public static double hdm_fromeXgetY(

       [ExcelArgument(Name = "地面数据", Description = "选择连续的三列数据 <桩号,偏距,高程>")] double[,] dmxy,
       [ExcelArgument(Name = "设计K", Description = "计算的设计桩号")] double k,
       [ExcelArgument(Name = "设计宽度", Description = "左负右正")] double b)
    {
        List<double> x = new List<double>();
        List<double> y = new List<double>();
        for (int i = 0; i < dmxy.GetLength(0); i++)
        {
            if (Math.Abs(k - dmxy[i,0]) <= 0.01)
            {
                x.Add(dmxy[i,1]);
                y.Add(dmxy[i,2]);
            }
        }
        return Math.Round( hdm.Interpolate(x.ToArray(), y.ToArray(), b),4);
    }
    #endregion

    #region start end
    public void AutoOpen()
    {
        MessageBox.Show(ExcelDnaUtil.XllPath + "");
        IntelliSenseServer.Install();
       
    }
    public void AutoClose()
    {
        IntelliSenseServer.Uninstall();
    }
    #endregion




}
