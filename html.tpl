{% extends 'full.tpl' %}

{% block header %}
{{ super() }}
<style>
div.prompt {display:none}
</style>

<script type="text/javascript">
    show=true;
    function toggle(){
        if (show){
            $('div.input').hide();
        }else{
            $('div.input').show();
        }
        show = !show
    }
    
    if(window.attachEvent) {
    window.attachEvent('onload', toggle);
    } else {
    if(window.onload) {
        var curronload = window.onload;
        var newonload = function() {
            curronload();
            toggle();
        };
        window.onload = newonload;
    } else {
        window.onload = toggle;
    }
}   
</script>
{% endblock header %}
