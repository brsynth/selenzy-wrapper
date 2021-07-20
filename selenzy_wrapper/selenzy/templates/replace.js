//for second row seq ID

var rows = document.getElementsByTagName('tr');

//initialize loop

for (var i = 1; i <= rows.length; i++)	{

	var ID = rows[i].getElementsByTagName('td')[0];

	var mylink = "//www.uniprot.org/uniprot/" + rows[i].getElementsByTagName('td')[0].innerHTML;

	var aTag = document.createElement('a');

	aTag.setAttribute('href', mylink);

	aTag.innerHTML = rows[i].getElementsByTagName('td')[0].innerHTML;

	ID.replaceChild(aTag, ID.childNodes[0]);

}
