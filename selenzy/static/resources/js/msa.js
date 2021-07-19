// Main jQuery call
$(document)
    .ready(
	function() {

	    $.ajax({	
		data : {
		    sessionid: JSON.stringify(sessionid),
		},
		type : 'POST',
		url : '/msa',
		success : function(serverdata) {
		    var data = JSON.parse(serverdata);
		    var seqs = msa.io.fasta.parse(data.msa);
		    var opts = {
			el: document.getElementById("msa"),
			seqs: seqs,
			vis: {
			    conserv: true,
			    overviewbox: false,
			    seqlogo: true,
			    metacell: false,
			    labelId: false,
			    labelCheckbox: false
			},
			zoomer : {
			    columnWidth: 15,
			    rowHeight: 20,
			    labelNameLength: 450,
			    labelPartLength: 200,
			    labelFontsize: "6px",
			    alignmentHeight: 1000
			},
			// smaller menu for JSBin
			menu: "small",
			bootstrapMenu: true
		    };
		    var m = msa(opts);
		    m.render();
		},
		error : function() {
	    console.log('Sorry, no luck');
		}
	    });
		
});
// To do: display tree 
