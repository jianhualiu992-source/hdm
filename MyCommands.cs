namespace hdm;
 
/**
    public static class MyCommands
    {
        // We make a command macro that can be run by:
     
        [ExcelCommand(MenuName = "hdm", MenuText = "Dump Data", ShortCut = "^D")]
        public static void DumpData()
        {
            // We always get the root Application object with a call to ExcelDnaUtil.Application
            // If we reference both Window Forms, and the Excel interop assemblies, 
            // we need to be a bit careful about which 'Application' we mean.
            Application app = (Application)ExcelDnaUtil.Application;
            var sheet = app.ActiveSheet;
            Range targetRange = sheet.Range["A1:C2"];
            object[,] newValues = new object[,] { { "One", 2, "Three" }, { true, System.DateTime.Now, "" } };
            targetRange.Value = newValues;

            // Apply some formatting, so that the time is displayed correctly
            Range dateCell = targetRange.Cells[2, 2];
            dateCell.NumberFormat = "hh:mm:ss";

        }
    }
  
**/
 
