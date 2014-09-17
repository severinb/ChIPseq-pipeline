$(function () {
    'use strict';

    jQuery.fn.dataTableExt.oSort['ratio-asc']  = function(x,y) {
        var a = eval(x);
        var b = eval(y);
        return ((a < b) ? -1 : ((a > b) ?  1 : 0));
    };
 
    jQuery.fn.dataTableExt.oSort['ratio-desc'] = function(x,y) {
        var a = eval(x);
        var b = eval(y);
        return ((a < b) ?  1 : ((a > b) ? -1 : 0));
    };

    /* peak table */
    var column_defs = [{"bSearchable": false, "sType": "numeric", "aTargets": [1]},
                       {"bSearchable": false, "sType": "numeric", "aTargets": [3]},
                       {"bSearchable": false, "sType": "numeric", "aTargets": [5]}]
    var columns = $("#peak_table th")
    var selector = "";
    var column_num = columns.length;
    for (var i = 6; i < column_num; i++) {
        column_defs.push({"bSearchable": false, "bVisible": false, "aTargets": [i]})
        selector += "<li data-column=\"" + i + "\" data-status=\"0\"><a>" + columns[i].textContent + " (show)</a></li>";
    }
    
    $("#motif_select").html(selector);
    $("#peak_table").dataTable({"sDom": "fliprt",
                                "iDisplayLength": 10,
                                   "aoColumnDefs": column_defs,
                                   "aaSorting": [ [1,'desc']]});
    
    $("#similarity_table").dataTable({"sDom": "lfiprt",
                                "iDisplayLength": 10,
                                "aoColumnDefs": [
                                    {"bSortable": false, "aTargets": [0]},
                                    {"bSearchable": false, "sType": "numeric", "aTargets": [1]},
                                    {"bSortable": false, "bSearchable": false, "aTargets": [2]}],
                                "aaSorting": [ [1,'asc']]});
    
     jQuery.fn.dataTableExt.oSort['peak-ratio-asc']  = function(x,y) {
        return ((eval(x) < eval(y)) ? -1 : ((eval(x) > eval(y)) ?  1 : 0));
     };
     jQuery.fn.dataTableExt.oSort['peak-ratio-desc']  = function(x,y) {
        return ((eval(x) < eval(y)) ? 1 : ((eval(x) > eval(y)) ?  -1 : 0));
     };

    $(".motif_table").dataTable({"sDom": "pirt",
                                "iDisplayLength": 10,
                                   "aoColumns": [
                                       {"bSortable": false},
                                       {"bSortable": false},
                                       {"bSortable": false},
                                       {"bSortable": false},
                                       {"bSortable": false},
                                       {"bSortable": false},
                                       {"bSortable": false},
                                       {"sType": "peak-ratio"}],
                                   "aaSorting": [ [2,'asc']]});

    $(".motif_tops_table").dataTable({"sDom": "pirt",
                                "iDisplayLength": 10,
                                   "aoColumns": [
                                       {"bSortable": false},
                                       {"bSortable": false},
                                       null,
                                       null,
                                       null,
                                       null,
                                       {"sType": "peak-ratio"}],
                                   "aaSorting": [ [2,'desc']]});
    

    $('ul.nav.nav-pills li a').click(function() {           
        $(this).parent().addClass('active').siblings().removeClass('active');           
    });
    $(".sh_switch").click(function(){
        var images = $(this).parent().parent().children(".images");
        if ($(images).css("display") == "block") {
            $(this).html("show");

        } else {
            $(this).html("hide");
        }
        $(images).toggle();
    });

    $("#motif_select li").click(function () {fhShowHide(this);});
});

function fhShowHide(el)
{
    /* Get the DataTables object again - this is not a recreation, just a get of the object */
    var oTable = $('#peak_table').dataTable();
    var column = $(el).data("column");

    /* get and change status */
    var status = $(el).data("status");
    console.log(status);
    $(el).data("status", Math.abs(status - 1));

    var text = $(el).children().html();
    if (status === 0) {
        text = text.replace(" (show)", " (hide)");
    }
    else {
        text = text.replace(" (hide)", " (show)");
    }
    $(el).children().html(text);
    
    var bVis = oTable.fnSettings().aoColumns[column].bVisible;
    oTable.fnSetColumnVis( column, bVis ? false : true );
}
