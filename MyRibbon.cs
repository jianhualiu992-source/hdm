
using hengduanmian;
namespace hdm
{
    [ComVisible(true)]
    public class MyRibbon : ExcelRibbon
    {
        public override string GetCustomUI(string RibbonID)
        {
            return RibbonResources.Ribbon;
        }

        public override object? LoadImage(string imageId)
        {
            // This will return the image resource with the name specified in the image='xxxx' tag
            return RibbonResources.ResourceManager.GetObject(imageId);
        }

        public void OnButtonPressed(IRibbonControl control)
        {
            string buttonID = control.Id;
            var acad = AcadHelper.GetActiveAutoCAD();
            acad.Visible = true;
            switch (buttonID)
            {
                #region 画圆
                case "addcircle":
                    double[] center = [];
                    double r = 0.2;
                    for (int i = 0; i < 100; i++)
                    {
                        center = [i, 0, 0];
                        acad.ActiveDocument.ModelSpace.AddCircle(center, r);
                    }

                    System.Windows.Forms.MessageBox.Show("ok!");
                    break;
                #endregion

                #region 获取点坐标
                case "getpoint":
                    //acad.ActiveDocument.ModelSpace.AddCircle(center, r);
                    Application app = (Application)ExcelDnaUtil.Application;
                    var cell = app.ActiveCell;
                    var hang = cell.Row;
                    var lie = cell.Column;
                    object basePoint = new double[] { 0, 0, 0 }; // 可以是 new double[] { 0, 0, 0 } 或 null
                    string prompt = "\n请选择一个点: ";
                    object result = acad.ActiveDocument.Utility.GetPoint(basePoint, prompt);
                    try
                    {
                        while (result != null)
                        {
                            double[] point = (double[])result;
                            app.Cells[hang, lie] = Math.Round(point[0], 3);
                            app.Cells[hang, lie + 1] = Math.Round(point[1], 3);
                            hang += 1;
                            result = acad.ActiveDocument.Utility.GetPoint(basePoint, prompt);
                        }
                    }
                    catch (Exception ex) { }
                    break;
                    #endregion




            }
        }
    }
}
