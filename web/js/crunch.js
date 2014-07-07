$(function () {
    'use strict';
    var files4upload = [];
    var uploaded_files = 0;

    $('#fileupload').bind('fileuploadadd', function (e, data) {
        console.log("add: "+ data.files[0].name)
        console.log("add: " + files4upload.indexOf(data.files[0].name))
        if (files4upload.indexOf(data.files[0].name) !== -1){
            data.files[0].error = "Duplicated";
        }
        console.log("add: "+ files4upload)
    });
    $('#fileupload').bind('fileuploadadded', function (e, data) {
        console.log("added: " + data.files[0].name)
        console.log("added: " + data.files[0].error)
        if (data.files[0].error === undefined) {
            files4upload.push(data.files[0].name);
        }
        console.log("added: "+ files4upload)
    });

    $('#fileupload').bind('fileuploadfailed', function (e, data) {
        if (data.files[0].error === null) {
            files4upload.splice(files4upload.indexOf(data.files[0].name), 1);
        }
    });
    
    $('#fileupload').bind('fileuploadstarted', function (e, data) {
        $(".cancel").hide();
    });

    $('#fileupload').bind('fileuploadcompleted', function (e, data) {
        uploaded_files += 1;
        console.log("uploaded files " + uploaded_files + " files4upload.length " + files4upload.length)
        if (uploaded_files == files4upload.length) {
            submit_user_data();
        }
    });

    $('#fileupload').fileupload({
        maxChunkSize: 10000000,
        acceptFileTypes: /(\.|\/)(fastq|gz)$/i
    });
    $("#optional").accordion({
	    collapsible: true,
	    heightStyle: "content"
		});

    $("#email").tooltip();
    $("#project").tooltip();
});



function submit_user_data() {
    var email = $("#email").val();
    var project = $("#project").val();
    var organism = $("#organism input:checked").val() === undefined ? "" : $("#organism input:checked").val();
    var fg_window = $("#fg_window").val();
    var bg_window = $("#bg_window").val();
    var step = $("#step").val();
    var fdr = $("#fdr").val();
    var adaptor = $("#adaptor").val();
    $.post("/fcgi/crunch.fcgi/run", {
        "email": email,
	"project": project,
        "organism": organism,
        "fg_window": fg_window,
        "bg_window": bg_window,
        "step": step,
	"fdr": fdr,
        "adaptor": adaptor
    },
           function(x){
               if(x.match(/^data_\S+$/) !== null) {
                   window.location = "http://" + window.location.hostname + "/CRUNCH/scratch/" + x + "/report";
               }
               else {
                   $(".container").html(x+"<br><br><h2><a href=''>Go back to main page</a></h2>");
               }
               $("#form").html("<br><br><h2>Done</h2>");
           });
}
