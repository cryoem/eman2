function listSearch() {
    inp = document.getElementById("input");
    filter = inp.value.toLowerCase();
    ul = document.getElementById("progsUL");
    li = ul.getElementsByTagName("li");
    for (let i = 0; i < li.length; ++i) {
        a = li[i].getElementsByTagName("a")[0];
        txt = a.innerText;
        if (txt.toLowerCase().indexOf(filter) > -1)
            li[i].style.display = "";
        else
            li[i].style.display = "none";
    }
}

document.addEventListener('keyup', function(event){
    if(event.key === "Escape") {
        document.getElementById("input").value="";
        listSearch()
    }
});
