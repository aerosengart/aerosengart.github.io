// Code from https://www.w3schools.com/howto/howto_js_snackbar.asp
function SnackbarFunc(id) {
    // Get the snackbar DIV
    var x = document.getElementById(id);
  
    // Add the "show" class to DIV
    x.className = "show";
  
    // After 3 seconds, remove the show class from DIV
    setTimeout(function(){ x.className = x.className.replace("show", ""); }, 3000);
}